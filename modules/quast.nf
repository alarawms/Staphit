nextflow.enable.dsl=2

process QUAST {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/quast", mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_quast_report"), path("${sample_id}_quast_metrics.tsv"), emit: quast_out

    script:
    """
    # Run QUAST with TSV output for easy parsing
    quast.py $assembly -o ${sample_id}_quast_report \
        --min-contig 1 \
        --output-format tsv \
        --threads ${task.cpus}

    # Extract key metrics to a simple TSV for filtering
    # QUAST creates a report.txt with all metrics
    REPORT="${sample_id}_quast_report/report.txt"
    
    # Extract metrics (QUAST report.txt format: Metric\\tValue)
    TOTAL_LENGTH=$(awk '/Total length/ {print $2}' "$REPORT")
    NUM_CONTIGS=$(awk '/# contigs/ {print $2}' "$REPORT")
    N50=$(awk '/N50/ {print $2}' "$REPORT")
    GC=$(awk '/GC %/ {print $2}' "$REPORT")
    
    # Calculate estimated completeness (genome fraction)
    # For S. aureus ~2.8Mb, but we can use genome fraction from QUAST
    GENOME_FRACTION=$(awk '/Genome fraction (%)/ {print $2}' "$REPORT")
    
    # Write metrics TSV
    echo -e "${sample_id}\\t${TOTAL_LENGTH}\\t${NUM_CONTIGS}\\t${N50}\\t${GC}\\t${GENOME_FRACTION}" > "${sample_id}_quast_metrics.tsv"
    """
}
