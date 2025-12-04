process SCCMEC {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/sccmec", mode: 'copy'
    container 'alarawms/sccmec_typer:latest'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("*.tsv"), emit: report
    tuple val(sample_id), path("*.json"), emit: json
    tuple val(sample_id), path("*.csv"), emit: csv
    path "versions.yml", emit: versions

    script:
    """
    python /app/bin/sccmec_typer.py \
        --1 ${assembly} \
        --db /app/db/sccmec_targets.fasta \
        --output ${sample_id}_sccmec \
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sccmec_typer: 1.0.0
    END_VERSIONS
    """
}
