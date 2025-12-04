nextflow.enable.dsl=2

process AMRFINDERPLUS {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/amrfinderplus", mode: 'copy'
    container 'staphb/ncbi-amrfinderplus:latest'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("*.tsv"), emit: report
    path "versions.yml", emit: versions

    script:
    """
    amrfinder -n ${assembly} --organism Staphylococcus_aureus --threads ${task.cpus} > ${sample_id}_amrfinder.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
    END_VERSIONS
    """
}
