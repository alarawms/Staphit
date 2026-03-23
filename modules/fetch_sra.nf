nextflow.enable.dsl=2

process FETCH_SRA {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/sra_downloads", mode: 'copy', pattern: "failed_sra.tsv"

    input:
    tuple val(sample_id), val(sra_id), val(reads_limit_param)

    output:
    tuple val(sample_id), path(fastq) emit: fastq
    path "failed_sra.tsv", emit: failures, optional: true

    script:
    def reads_limit_cmd = reads_limit_param ? "-X ${reads_limit_param}" : ""
    def max_retries = params.max_sra_retries ?: 3
    
    """
    # Configure sra-tools
    mkdir -p \$HOME/.ncbi
    if [ -f /root/.ncbi/user-settings.mkfg ]; then
        cp /root/.ncbi/user-settings.mkfg \$HOME/.ncbi/user-settings.mkfg
    else
        printf '/LIBS/GUID = "88608c0e-f078-4908-bf3d-94483989301a"\\n' > \$HOME/.ncbi/user-settings.mkfg
        printf '/libs/cloud/report_instance_identity = "true"\\n' >> \$HOME/.ncbi/user-settings.mkfg
    fi

    mkdir fastq
    
    # Prefetch with retry (exponential backoff)
    n=0
    PREFETCH_SUCCESS=false
    while [ "\$n" -lt ${max_retries} ]; do
        if prefetch "$sra_id" 2>&1; then
            PREFETCH_SUCCESS=true
            break
        else
            n=\$((n+1))
            WAIT_TIME=\$((5 * 2 ** (n-1)))
            echo "Prefetch failed (attempt \$n/${max_retries}), retrying in ${WAIT_TIME}s..."
            sleep $WAIT_TIME
        fi
    done
    
    if [ "\$PREFETCH_SUCCESS" = false ]; then
        echo "ERROR: Prefetch failed after ${max_retries} retries for $sra_id"
        echo -e "${sample_id}\\t$sra_id\\tPREFETCH_FAILED\\tExhausted ${max_retries} retries" >> failed_sra.tsv
        exit 1
    fi

    # Use fastq-dump with error handling
    if fastq-dump --split-3 --gzip $reads_limit_cmd --outdir fastq $sra_id 2>&1; then
        echo "fastq-dump succeeded for $sra_id"
    else
        echo "ERROR: fastq-dump failed for $sra_id"
        echo -e "${sample_id}\\t$sra_id\\tFASTQ_DUMP_FAILED\\tSee error above" >> failed_sra.tsv
        exit 1
    fi
    
    # Verify paired-end output
    if [ -f fastq/${sra_id}_1.fastq.gz ] && [ -f fastq/${sra_id}_2.fastq.gz ]; then
        mv fastq/${sra_id}_1.fastq.gz fastq/${sample_id}_1.fastq.gz
        mv fastq/${sra_id}_2.fastq.gz fastq/${sample_id}_2.fastq.gz
    else
        echo "ERROR: Expected paired-end files (_1.fastq.gz and _2.fastq.gz) not found"
        echo -e "${sample_id}\\t$sra_id\\tPAIRED_END_MISSING\\tfastq-dump did not produce PE output" >> failed_sra.tsv
        ls -R fastq/ >> failed_sra.tsv
        exit 1
    fi
    """
}
