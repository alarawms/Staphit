nextflow.enable.dsl=2

process MLST {
    tag "$sample_id"
    publishDir "${params.outdir}/mlst/$sample_id", mode: 'copy'

    input:
    tuple val(sample_id), path(scaffolds)

    output:
    tuple val(sample_id), path("*.tsv")

    script:
    """
    mlst --csv $scaffolds > ${sample_id}.tsv
    """
}
