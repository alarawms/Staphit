nextflow.enable.dsl=2

process FETCH_SRA {
    tag "$sample_id"
    label 'process_low'

    input:
    tuple val(sample_id), val(sra_id), val(reads_limit_param)

    output:
    tuple val(sample_id), path(fastq)

    script:
    def reads_limit_cmd = reads_limit_param ? "-X ${reads_limit_param}" : ""
    """
    # Configure sra-tools
    mkdir -p \$HOME/.ncbi
    if [ -f /root/.ncbi/user-settings.mkfg ]; then
        cp /root/.ncbi/user-settings.mkfg \$HOME/.ncbi/user-settings.mkfg
    else
        # Create a minimal config with a placeholder GUID
        printf '/LIBS/GUID = "88608c0e-f078-4908-bf3d-94483989301a"\\n' > \$HOME/.ncbi/user-settings.mkfg
        printf '/libs/cloud/report_instance_identity = "true"\\n' >> \$HOME/.ncbi/user-settings.mkfg
    fi

    mkdir fastq
    
    # Retry prefetch a few times if it fails
    n=0
    until [ "\$n" -ge 3 ]
    do
       prefetch $sra_id && break
       n=\$((n+1)) 
       sleep 5
       echo "Prefetch failed, retrying (\$n/3)..."
    done

    # Use fastq-dump which is more reliable for splitting and supports gzip directly
    # Limit spots based on params.reads_limit if set
    fastq-dump --split-3 --gzip $reads_limit_cmd --outdir fastq $sra_id
    
    echo "Listing fastq directory content:"
    ls -l fastq/

    # Rename to match expected output
    # Check if paired files exist, otherwise fail (pipeline expects PE)
    if [ -f fastq/${sra_id}_1.fastq.gz ]; then
        mv fastq/${sra_id}_1.fastq.gz fastq/${sample_id}_1.fastq.gz
        mv fastq/${sra_id}_2.fastq.gz fastq/${sample_id}_2.fastq.gz
    else
        echo "Error: Expected paired-end output (_1.fastq.gz and _2.fastq.gz) but not found."
        # List what was found for debugging
        ls -R fastq/
        exit 1
    fi
    """
}
