process AGR_TYPING {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/agr_typing", mode: 'copy'
    containerOptions '--entrypoint ""'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("agr_results/"), emit: results
    tuple val(sample_id), path("agr_results/result.json"), emit: report

    script:
    """
    staph_agr_typer run --fasta ${assembly} -o agr_results || true

    # Ensure result.json exists even if typing failed
    if [ ! -f agr_results/result.json ]; then
        mkdir -p agr_results
        echo '{"agr_group": "ND", "confidence": 0.0}' > agr_results/result.json
    fi
    """
}
