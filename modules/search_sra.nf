nextflow.enable.dsl=2

process SEARCH_SRA {
    container 'ncbi/edirect:latest'
    publishDir '.', mode: 'copy', overwrite: true

    input:
    val species
    val limit

    output:
    path "samplesheet.csv", emit: samplesheet

    shell:
    '''
    # Search for paired-end Illumina runs for the species
    query="!{species} AND paired[Layout] AND illumina[Platform]"
    
    echo "Searching NCBI SRA for: $query"
    
    # 1. Fetch RunInfo to find valid runs (spots > 0)
    # Fetch 1000 to ensure we find valid runs
    esearch -db sra -query "$query" | \
    efetch -format runinfo -stop 1000 > runinfo.csv

    if [ ! -s runinfo.csv ]; then
        echo "No results found for $query"
        exit 1
    fi

    # 2. Filter for valid runs and select top N
    # Extract Run IDs of valid runs (spots > 0)
    tail -n +2 runinfo.csv | awk -F ',' '$4 > 0 {print $1}' | head -n !{limit} > selected_runs.txt
    
    count=$(wc -l < selected_runs.txt)
    if [ "$count" -eq 0 ]; then
        echo "No valid runs (spots > 0) found."
        exit 1
    fi
    echo "Selected $count runs for processing."

    # 3. Generate samplesheet.csv
    echo "sample,sra,fastq_1,fastq_2" > samplesheet.csv
    while read run; do
        echo "$run,$run,," >> samplesheet.csv
    done < selected_runs.txt

    echo "Generated samplesheet.csv"
    '''
}
