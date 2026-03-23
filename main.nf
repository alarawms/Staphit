#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// -- PARAMETERS --
params.outdir = './results'
params.input = 'samplesheet.csv'

// -- MODULES --
include { FETCH_SRA } from './modules/fetch_sra.nf'
include { TRIMMOMATIC } from './modules/trimmomatic.nf'
include { FASTQC } from './modules/fastqc.nf'
include { SPADES } from './modules/spades.nf'
include { SKESA } from './modules/skesa.nf'
include { PROKKA } from './modules/prokka.nf'
include { ABRICATE } from './modules/abricate.nf'
include { MLST } from './modules/mlst.nf'
include { MULTIQC } from './modules/multiqc.nf'
include { PANAROO } from './modules/panaroo.nf'
include { IQTREE } from './modules/iqtree.nf'
include { SNP_DISTS } from './modules/snp_dists.nf'
include { QUAST } from './modules/quast.nf'
include { MASH } from './modules/mash.nf'
include { AMRFINDERPLUS } from './modules/amrfinderplus.nf'
include { SPATYPER } from './modules/spatyper.nf'
include { SCCMEC } from './modules/sccmec.nf'
include { SCREEN_SPECIES } from './modules/screen_species.nf'
include { QC_FILTER } from './modules/qc_filter.nf'
include { GENERATE_SUMMARY } from './modules/generate_summary.nf'
include { AGR_TYPING } from './modules/agr_typing.nf'
include { FETCH_RESFINDER_DB; INDEX_DB; KMA } from './modules/kma.nf'
include { VISUALIZATION } from './modules/visualization.nf'

include { SEARCH_SRA } from './modules/search_sra.nf'
include { FETCH_METADATA } from './modules/fetch_metadata.nf'
include { AGGREGATOR } from './modules/aggregator.nf'
include { SUMMARY_MERGER } from './modules/summary_merger.nf'

// -- WORKFLOW --
workflow {
    main:
        // Debug prints
        log.info "Params: species=${params.species}, search_limit=${params.search_limit}, input=${params.input}, assembler=${params.assembler}"

        // Initialize metadata channel
        ch_metadata_json = Channel.empty()

        if (params.species) {
            log.info "Running SRA Search for species: ${params.species}"
            SEARCH_SRA(params.species, params.search_limit)
            input_csv = SEARCH_SRA.out.samplesheet
            
            // Fetch rich metadata using Python module
            FETCH_METADATA(input_csv)
            ch_metadata_json = FETCH_METADATA.out.metadata_json
        } else {
            log.info "Using provided samplesheet: ${params.input}"
            // If not searching, we expect the file to exist.
            input_csv = Channel.fromPath(params.input, checkIfExists: true)
            
            // Try to find metadata.json if available locally, otherwise create empty one
            if (file('metadata.json').exists()) {
                ch_metadata_json = Channel.fromPath('metadata.json')
            } else {
                // Create a dummy metadata file so the aggregator can still run
                def dummy = file("${workDir}/metadata.json")
                dummy.text = '[]'
                ch_metadata_json = Channel.of(dummy)
            }
        }

        // Create a channel from the sample sheet
        ch_input = input_csv
            .splitCsv(header: true)
            .map { row ->
                def meta = [id: row.sample]
                if (row.sra && !row.sra.trim().isEmpty()) {
                    return [ meta, row.sra, null, null ] // SRA sample
                } else if (row.fastq_1 && !row.fastq_1.trim().isEmpty() && row.fastq_2 && !row.fastq_2.trim().isEmpty()) {
                    return [ meta, null, file(row.fastq_1), file(row.fastq_2) ] // Local sample
                } else {
                    error "Invalid samplesheet entry or incomplete FastQ paths: ${row}"
                }
            }
            .branch {
                sra_samples: it[1] != null
                local_samples: it[2] != null
            }
            .set { ch_processed_samples }

        // Apply max_downloads limit to SRA samples if set
        ch_sra_to_fetch = ch_processed_samples.sra_samples
        if (params.max_downloads != null && params.max_downloads as int > 0) {
            ch_sra_to_fetch = ch_sra_to_fetch | take(params.max_downloads as int)
        }

        // --- SRA Processing ---
        FETCH_SRA(ch_sra_to_fetch.map { row_tuple -> [row_tuple[0].id, row_tuple[1], params.reads_limit] })
        
        reads_from_sra_ch = FETCH_SRA.out.map { sample_id, path ->
            def read1 = path.resolve(sample_id + "_1.fastq.gz")
            def read2 = path.resolve(sample_id + "_2.fastq.gz")
            [ sample_id, [read1, read2] ]
        }

        // --- Local Files Processing ---
        reads_from_local_ch = ch_processed_samples.local_samples.map { meta, _, fq1, fq2 ->
            [ meta.id, [fq1, fq2] ]
        }

        // --- Mix SRA and Local Reads ---
        reads_ch = reads_from_sra_ch.mix(reads_from_local_ch)

        // --- QC with Trimmomatic and FastQC ---
        ch_adapters = Channel.fromPath("${projectDir}/assets/TruSeq3-PE.fa", checkIfExists: true)
        TRIMMOMATIC(reads_ch, ch_adapters)
        FASTQC(TRIMMOMATIC.out.trimmed_reads)
        
        ch_trimmed_reads = TRIMMOMATIC.out.trimmed_reads.map { id, files -> [id, files] }
        
        // --- Assembly ---
        if (params.assembler == 'skesa') {
            SKESA(ch_trimmed_reads)
            ch_raw_assemblies = SKESA.out
        } else {
            SPADES(ch_trimmed_reads)
            ch_raw_assemblies = SPADES.out
        }

        // Filter out poor assemblies (< 500 KB = junk for S. aureus ~2.8 Mb)
        ch_assemblies = ch_raw_assemblies.filter { sample_id, fasta ->
            fasta.size() > params.min_assembly_size
        }

        // --- Species Screening ---
        if (!params.skip_species_check) {
            log.info "Running species screening against S. aureus USA300 (mash dist < ${params.min_mash_dist})..."
            SCREEN_SPECIES(ch_assemblies)
            ch_screened_assemblies = SCREEN_SPECIES.out.screened_assemblies
                .filter { sample_id, fasta, screen_file ->
                    def pass = new File(screen_file).text.split('\t')[1] == 'PASS'
                    if (!pass) {
                        log.warn "Sample ${sample_id} failed species screen - filtered out"
                    }
                    return pass
                }
                .map { sample_id, fasta, screen_file -> [sample_id, fasta] }
            
            # Copy screen summary to output
            SCREEN_SPECIES.out.summary.map { it -> file(it).copyTo("${params.outdir}/screen_species_summary.tsv") }
        } else {
            log.info "Skipping species check (--skip_species_check set)"
            ch_screened_assemblies = ch_assemblies
        }

        // --- QC with QUAST and Filtering ---
        QUAST(ch_screened_assemblies)
        
        // Apply QC filters if any are set
        if (params.min_completeness || params.max_contamination || params.min_n50 || params.max_n50) {
            log.info "Applying QC filters: completeness=${params.min_completeness}%, contamination=${params.max_contamination}%, N50=${params.min_n50}-${params.max_n50}"
            QC_FILTER(QUAST.out.quast_out)
            ch_qc_filtered = QC_FILTER.out.qc_filtered
                .filter { sample_id, quast_dir, qc_file ->
                    def pass = new File(qc_file).text.split('\t')[1] == 'PASS'
                    if (!pass) {
                        reason = new File(qc_file).text.split('\t')[8]
                        log.warn "Sample ${sample_id} failed QC: ${reason}"
                    }
                    return pass
                }
            # For aggregation, we need: sample_id -> quast_dir (keeps metrics inside)
            ch_quast_for_agg = QUAST.out.quast_out.map { sample_id, quast_dir, metrics -> [sample_id, quast_dir] }
        } else {
            log.info "No QC filters applied (all QC params are null)"
            ch_qc_filtered = QUAST.out.quast_out.map { sample_id, quast_dir, metrics_file -> [sample_id, quast_dir] }
            ch_quast_for_agg = QUAST.out.quast_out.map { sample_id, quast_dir, metrics_file -> [sample_id, quast_dir] }
        }

        MASH(ch_qc_filtered)
        AMRFINDERPLUS(ch_qc_filtered)
        SPATYPER(ch_qc_filtered)
        SCCMEC(ch_qc_filtered)
        AGR_TYPING(ch_qc_filtered)
        PROKKA(ch_qc_filtered)
        ABRICATE(ch_qc_filtered)
        MLST(ch_qc_filtered)

        // --- Read-based AMR (KMA) ---
        FETCH_RESFINDER_DB()
        INDEX_DB(FETCH_RESFINDER_DB.out.db)
        KMA(ch_trimmed_reads, INDEX_DB.out.indexed_db.collect())

        // --- Aggregation ---
        // Prepare inputs:
        // Trimmomatic log and FastQC output
        ch_trim_log = TRIMMOMATIC.out.log
        ch_fastqc_out = FASTQC.out.map { id, files -> [id, files] }
        
        // Join all outputs by sample_id
        ch_agg_in = ch_trim_log
            .join(ch_fastqc_out)
            .join(ch_quast_for_agg)
            .join(MLST.out)
            .join(ABRICATE.out)
            .join(AMRFINDERPLUS.out.report)
            .join(MASH.out.sketch)
            .join(SPATYPER.out.report)
            .join(SCCMEC.out.report)
            .join(AGR_TYPING.out.report)
            .join(KMA.out.results)
            
        // Combine with metadata.json
        ch_agg_final = ch_agg_in.combine(ch_metadata_json)
        
        AGGREGATOR(ch_agg_final)
        
        // Collect all summary CSVs and merge them
        SUMMARY_MERGER(AGGREGATOR.out[1].collect())

        // --- Visualization ---
        VISUALIZATION(SUMMARY_MERGER.out)

        // Pangenome and Phylogeny
        PANAROO(PROKKA.out.collect())
        IQTREE(PANAROO.out.aln)
        SNP_DISTS(PANAROO.out.aln)

        // Collect all the outputs and pass them to MultiQC
        ch_multiqc_in = channel.empty()
        ch_multiqc_in = ch_multiqc_in.mix(FASTQC.out.map{ it[1] }.collect())
        ch_multiqc_in = ch_multiqc_in.mix(TRIMMOMATIC.out.log.map{ it[1] }.collect())
        ch_multiqc_in = ch_multiqc_in.mix(ch_assemblies.map{ it[1] }.collect())
        ch_multiqc_in = ch_multiqc_in.mix(PROKKA.out.collect())
        ch_multiqc_in = ch_multiqc_in.mix(ABRICATE.out.map{ it[1] }.collect())
        ch_multiqc_in = ch_multiqc_in.mix(MLST.out.map{ it[1] }.collect())
        ch_multiqc_in = ch_multiqc_in.mix(ch_qc_filtered.map{ it[1] }.collect())
        
        # Add QC filter summary if filters were applied
        if (params.min_completeness || params.max_contamination || params.min_n50 || params.max_n50) {
            QC_FILTER.out.summary.map { it -> file(it).copyTo("${params.outdir}/qc_summary.tsv") }
        }
        
        MULTIQC(ch_multiqc_in.collect())

        // --- HTML Summary Dashboard ---
        // Prepare optional summaries
        ch_screen_summary = Channel.empty()
        ch_qc_summary = Channel.empty()
        
        if (!params.skip_species_check) {
            ch_screen_summary = SCREEN_SPECIES.out.summary
        }
        
        if (params.min_completeness || params.max_contamination || params.min_n50 || params.max_n50) {
            ch_qc_summary = QC_FILTER.out.summary
        }
        
        GENERATE_SUMMARY(
            ch_screen_summary.collect(),
            ch_qc_summary.collect(),
            ch_quast_for_agg.map { it[1] }.collect(),
            MLST.out,
            SPATYPER.out.report,
            SCCMEC.out.report,
            AGR_TYPING.out.report,
            ABRICATE.out,
            AMRFINDERPLUS.out.report
        )
}

