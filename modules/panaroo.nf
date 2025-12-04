nextflow.enable.dsl=2

process PANAROO {
    publishDir "${params.outdir}/panaroo", mode: 'copy'

    input:
    path gffs

    output:
    path "core_gene_alignment.aln", optional: true, emit: aln
    path "gene_presence_absence.csv", optional: true
    path "pan_genome_reference.fa", optional: true

    script:
    """
    # Panaroo requires at least 2 samples. If we have fewer, skip or mock.
    count=\$(ls *.gff | wc -l)
    if [ "\$count" -lt 2 ]; then
        echo "Not enough samples for Panaroo (needs > 1). Skipping."
        exit 0
    fi

    panaroo -i *.gff -o . --clean-mode strict --remove-invalid-genes -a core --aligner mafft -t ${task.cpus}
    """
}
