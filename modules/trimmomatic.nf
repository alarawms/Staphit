nextflow.enable.dsl=2

process TRIMMOMATIC {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.trimmed.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("*.trim.log"), emit: log

    script:
    def (r1, r2) = reads
    """
    trimmomatic PE -phred33 \\
        $r1 $r2 \\
        ${sample_id}_1.trimmed.fastq.gz ${sample_id}_1.unpaired.fastq.gz \\
        ${sample_id}_2.trimmed.fastq.gz ${sample_id}_2.unpaired.fastq.gz \\
        ILLUMINACLIP:$baseDir/assets/TruSeq3-PE.fa:2:30:10 \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> ${sample_id}.trim.log
    """
}

