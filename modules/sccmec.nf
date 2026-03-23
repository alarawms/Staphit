process SCCMEC {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/sccmec", mode: 'copy'
    container 'alarawms/sccmec_typer:latest'
    containerOptions '--entrypoint ""'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("*.tsv"), emit: report
    tuple val(sample_id), path("*.json"), emit: json
    tuple val(sample_id), path("*.csv"), emit: csv
    path "versions.yml", emit: versions

    script:
    """
    conda run --no-capture-output -n sccmec_typer \
        python /app/bin/sccmec_typer.py \
        --1 ${assembly} \
        -d /app/db/sccmec_targets.fasta \
        -o ${sample_id}_sccmec \
        --threads ${task.cpus} \
        --no-viz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sccmec_typer: 1.0.0
    END_VERSIONS
    """
}
