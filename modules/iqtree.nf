nextflow.enable.dsl=2

process IQTREE {
    publishDir "${params.outdir}/iqtree", mode: 'copy'

    input:
    path alignment

    output:
    path "*.treefile", optional: true
    path "*.iqtree", optional: true

    script:
    """
    # Check number of sequences in alignment (count '>' lines)
    seq_count=\$(grep -c "^>" $alignment)
    if [ "\$seq_count" -lt 3 ]; then
        echo "Alignment has only \$seq_count sequences. IQ-TREE requires at least 3. Skipping."
        exit 0
    fi

    iqtree -s $alignment -m MFP -nt AUTO -ntmax ${task.cpus} -bb 1000
    """
}
