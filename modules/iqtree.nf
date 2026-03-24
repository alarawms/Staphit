nextflow.enable.dsl=2

process IQTREE {
    publishDir "${params.outdir}/iqtree", mode: 'copy'

    input:
    path alignment
    path seed_tree, stageAs: 'seed.treefile'

    output:
    path "*.treefile", optional: true
    path "*.iqtree", optional: true

    script:
    def seed_flag = seed_tree.name != 'NO_SEED_TREE' ? "-t seed.treefile" : ''
    """
    # Check number of sequences in alignment (count '>' lines)
    seq_count=\$(grep -c "^>" $alignment)
    if [ "\$seq_count" -lt 3 ]; then
        echo "Alignment has only \$seq_count sequences. IQ-TREE requires at least 3. Skipping."
        exit 0
    fi

    iqtree2 -s $alignment \
        -m ${params.iqtree_model ?: 'GTR+F+I'} \
        -nt AUTO -ntmax ${task.cpus} \
        -bb ${params.iqtree_bb ?: 1000} \
        ${params.iqtree_fast ? '-fast' : ''} \
        ${seed_flag}
    """
}
