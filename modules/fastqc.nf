nextflow.enable.dsl=2

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("fastqc_reports/*")

    script:
    """
    mkdir fastqc_reports
    fastqc -o fastqc_reports $reads
    """
}
