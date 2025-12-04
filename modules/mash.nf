nextflow.enable.dsl=2

process MASH {
    tag "$sample_id"
    label 'process_medium'
    container 'staphb/mash:latest'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("*.msh"), emit: sketch
    tuple val(sample_id), path("*.dist"), emit: dist
    path "versions.yml", emit: versions

    script:
    """
    mash sketch -p ${task.cpus} -o ${sample_id} ${assembly}
    
    # Calculate distance against itself (as a placeholder/sanity check) 
    # In a real scenario, you'd compare against a reference database
    mash dist ${sample_id}.msh ${sample_id}.msh > ${sample_id}.dist

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version)
    END_VERSIONS
    """
}
