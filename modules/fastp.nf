process FASTP {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_1.trim.fastq.gz"), path("${sample_id}_2.trim.fastq.gz"), emit: reads
    tuple val(sample_id), path("*.json"), emit: json
    tuple val(sample_id), path("*.html"), emit: html

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${sample_id}_1.trim.fastq.gz \\
        --out2 ${sample_id}_2.trim.fastq.gz \\
        --json ${sample_id}.fastp.json \\
        --html ${sample_id}.fastp.html \\
        --detect_adapter_for_pe \\
        --thread ${task.cpus}
    """
}
