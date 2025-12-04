nextflow.enable.dsl=2

process PROKKA {
    tag "$sample_id"
    publishDir "${params.outdir}/prokka/$sample_id", mode: 'copy'

    input:
    tuple val(sample_id), path(scaffolds)

    output:
    path("*.gff")

    script:
    """
    prokka --outdir . --force --prefix $sample_id --kingdom Bacteria --genus Staphylococcus --species aureus --strain $sample_id $scaffolds
    """
}
