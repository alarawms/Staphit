#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// -- PARAMETERS --
params.outdir = './results'
params.input = 'samplesheet.csv'

// -- MODULES --
include { FETCH_SRA } from './modules/fetch_sra.nf'
include { FASTP } from './modules/fastp.nf'
include { SPADES } from './modules/spades.nf'
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
        log.info "Params: species=${params.species}, search_limit=${params.search_limit}, input=${params.input}"

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
            
            // Try to find metadata.json if available locally
            if (file('metadata.json').exists()) {
                ch_metadata_json = Channel.fromPath('metadata.json')
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

        // --- QC with fastp ---
        FASTP(reads_ch)
        
        ch_trimmed_reads = FASTP.out.reads.map { id, fq1, fq2 -> [id, [fq1, fq2]] }
        
        // --- Assembly ---
        SPADES(ch_trimmed_reads)
        QUAST(SPADES.out)
        MASH(SPADES.out)
        AMRFINDERPLUS(SPADES.out)
        SPATYPER(SPADES.out)
        SCCMEC(SPADES.out)
        AGR_TYPING(SPADES.out)
        PROKKA(SPADES.out)
        ABRICATE(SPADES.out)
        MLST(SPADES.out)

        // --- Read-based AMR (KMA) ---
        FETCH_RESFINDER_DB()
        INDEX_DB(FETCH_RESFINDER_DB.out.db)
        KMA(ch_trimmed_reads, INDEX_DB.out.indexed_db.collect())

        // --- Aggregation ---
        // Prepare inputs:
        // FASTP output: [id, json, html] -> [id, json]
        ch_fastp_json = FASTP.out.json.map { id, json -> [id, json] }
        
        // Join all outputs by sample_id
        ch_agg_in = ch_fastp_json
            .join(QUAST.out)
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
        ch_multiqc_in = ch_multiqc_in.mix(FASTP.out.json.map{ it[1] }.collect())
        ch_multiqc_in = ch_multiqc_in.mix(SPADES.out.map{ it[1] }.collect())
        ch_multiqc_in = ch_multiqc_in.mix(PROKKA.out.collect())
        ch_multiqc_in = ch_multiqc_in.mix(ABRICATE.out.map{ it[1] }.collect())
        ch_multiqc_in = ch_multiqc_in.mix(MLST.out.map{ it[1] }.collect())
        ch_multiqc_in = ch_multiqc_in.mix(QUAST.out.map{ it[1] }.collect())
        ch_multiqc_in = ch_multiqc_in.mix(MASH.out.sketch.map{ it[1] }.collect())
        
        MULTIQC(ch_multiqc_in.collect())
}

