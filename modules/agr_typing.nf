process AGR_TYPING {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/agr_typing", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("agrvate_results/"), emit: results
    tuple val(sample_id), path("agrvate_results/result.json"), emit: report

    script:
    """
    staph_agr_typer run --fasta ${assembly} -o agrvate_results
    """
}
