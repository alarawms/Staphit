process FETCH_RESFINDER_DB {
    label 'process_low'
    
    output:
    path("resfinder_db"), emit: db

    script:
    """
    git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git
    """
}

process INDEX_DB {
    label 'process_medium'

    input:
    path(db)

    output:
    path("resfinder_db_indexed"), emit: indexed_db

    script:
    """
    mkdir resfinder_db_indexed
    cp -r ${db}/* resfinder_db_indexed/
    cd resfinder_db_indexed
    kma_index -i *.fsa -o resfinder_kma
    """
}

process KMA {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/kma", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path(db_dir)

    output:
    tuple val(sample_id), path("${sample_id}.res"), emit: results
    tuple val(sample_id), path("${sample_id}.mapstat"), emit: mapstat

    script:
    """
    kma -ipe ${reads[0]} ${reads[1]} -t_db ${db_dir}/resfinder_kma -o ${sample_id} -1t1
    [ -f "${sample_id}.res" ] || touch "${sample_id}.res"
    [ -f "${sample_id}.mapstat" ] || touch "${sample_id}.mapstat"
    """
}
