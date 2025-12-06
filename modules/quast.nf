nextflow.enable.dsl=2

process QUAST {
    tag "$sample_id"
    publishDir "${params.outdir}/quast", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_quast_report")

    script:
    """
    quast.py $assembly -o ${sample_id}_quast_report --min-contig 1
    """
}
