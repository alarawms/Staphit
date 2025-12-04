nextflow.enable.dsl=2

process ABRICATE {
    tag "$sample_id"
    publishDir "${params.outdir}/abricate/$sample_id", mode: 'copy'

    input:
    tuple val(sample_id), path(scaffolds)

    output:
    tuple val(sample_id), path("*.tab")

    script:
    """
    abricate --db resfinder $scaffolds > ${sample_id}_resfinder.tab
    abricate --db vfdb $scaffolds > ${sample_id}_vfdb.tab
    abricate --db plasmidfinder $scaffolds > ${sample_id}_plasmidfinder.tab
    """
}
