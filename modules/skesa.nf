nextflow.enable.dsl=2

process SKESA {
    tag "$sample_id"
    publishDir "${params.outdir}/skesa", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.scaffolds.fasta")

    script:
    def (r1, r2) = reads
    """
    skesa --reads $r1,$r2 --cores ${task.cpus} --memory ${task.memory.giga.intValue()} --contigs_out ${sample_id}.scaffolds.fasta
    """
}
