nextflow.enable.dsl=2

process SPATYPER {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/spatyper", mode: 'copy'
    container 'python:3.9-slim'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("*.tsv"), emit: report
    path "versions.yml", emit: versions

    script:
    """
    pip install --quiet spaTyper
    spaTyper -f ${assembly} --output ${sample_id}_spa.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spatyper: \$(spaTyper --version 2>&1 || echo "unknown")
    END_VERSIONS
    """
}
