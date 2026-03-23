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
    # Remove poor assemblies — S. aureus should have ~2500 CDS; anything under 500 is junk
    for f in *.gff; do
        cds_count=\$(grep -c "CDS" "\$f" 2>/dev/null || echo 0)
        if [ ! -s "\$f" ] || [ "\$cds_count" -lt 500 ]; then
            echo "Removing invalid GFF (\$cds_count CDS): \$f"
            rm -f "\$f"
        fi
    done

    count=\$(ls *.gff 2>/dev/null | wc -l)
    if [ "\$count" -lt 2 ]; then
        echo "Not enough valid samples for Panaroo (needs >= 2, found \$count). Skipping."
        exit 0
    fi

    panaroo -i *.gff -o . --clean-mode strict --remove-invalid-genes -a core --aligner mafft -t ${task.cpus}
    """
}
