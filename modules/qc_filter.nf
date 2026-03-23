nextflow.enable.dsl=2

process QC_FILTER {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/qc_filter", mode: 'copy'

    input:
    tuple val(sample_id), path(quast_dir), path(metrics_tsv)

    output:
    tuple val(sample_id), path(quast_dir), path("${sample_id}_qc.tsv"), emit: qc_filtered
    path "qc_summary.tsv", emit: summary

    script:
    """
    # Read metrics
    TOTAL_LENGTH=$(cut -f2 "$metrics_tsv")
    NUM_CONTIGS=$(cut -f3 "$metrics_tsv")
    N50=$(cut -f4 "$metrics_tsv")
    GC=$(cut -f5 "$metrics_tsv")
    GENOME_FRACTION=$(cut -f6 "$metrics_tsv")
    
    # Apply QC thresholds (only if set)
    QC_PASS=true
    REASON=""
    
    # Min completeness (via genome fraction)
    if [ -n "${params.min_completeness}" ]; then
        if (( $(echo "$GENOME_FRACTION < ${params.min_completeness}" | bc -l) )); then
            QC_PASS=false
            REASON="${REASON}Low_completeness(${GENOME_FRACTION}%) "
        fi
    fi
    
    # Max contamination (1 - genome_fraction)
    if [ -n "${params.max_contamination}" ]; then
        CONTAM=$(echo "100 - $GENOME_FRACTION" | bc -l)
        if (( $(echo "$CONTAM > ${params.max_contamination}" | bc -l) )); then
            QC_PASS=false
            REASON="${REASON}High_contamination(${CONTAM}%) "
        fi
    fi
    
    # Min N50
    if [ -n "${params.min_n50}" ]; then
        if (( $(echo "$N50 < ${params.min_n50}" | bc -l) )); then
            QC_PASS=false
            REASON="${REASON}Low_N50(${N50}) "
        fi
    fi
    
    # Max N50
    if [ -n "${params.max_n50}" ]; then
        if (( $(echo "$N50 > ${params.max_n50}" | bc -l) )); then
            QC_PASS=false
            REASON="${REASON}High_N50(${N50}) "
        fi
    fi
    
    # Clean reason
    REASON=$(echo "$REASON" | sed 's/ *$//')
    
    # Write QC result
    if [ "$QC_PASS" = true ]; then
        echo -e "${sample_id}\\tPASS\\t${TOTAL_LENGTH}\\t${NUM_CONTIGS}\\t${N50}\\t${GC}\\t${GENOME_FRACTION}\\t-" > "${sample_id}_qc.tsv"
    else
        echo -e "${sample_id}\\tFAIL\\t${TOTAL_LENGTH}\\t${NUM_CONTIGS}\\t${N50}\\t${GC}\\t${GENOME_FRACTION}\\t${REASON}" > "${sample_id}_qc.tsv"
    fi
    
    # Append to summary
    echo -e "${sample_id}\\t${QC_PASS}\\t${TOTAL_LENGTH}\\t${NUM_CONTIGS}\\t${N50}\\t${GC}\\t${GENOME_FRACTION}\\t${REASON}" >> qc_summary.tsv
    """
}
