nextflow.enable.dsl=2

process SPADES {
    tag "$sample_id"
    publishDir "${params.outdir}/spades", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.scaffolds.fasta")

    script:
    def (r1, r2) = reads
    """
    spades.py --only-assembler -1 $r1 -2 $r2 -o . -m ${task.memory.giga.intValue()}
    mv scaffolds.fasta ${sample_id}.scaffolds.fasta
    """
}
