nextflow.enable.dsl=2

process SNP_DISTS {
    publishDir "${params.outdir}/snp_dists", mode: 'copy'

    input:
    path alignment

    output:
    path "snp_distances.tsv"

    script:
    """
    snp-dists $alignment > snp_distances.tsv
    """
}
