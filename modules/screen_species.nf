process SCREEN_SPECIES {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/screen_species", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path(fasta), path("${sample_id}.screen.tsv"), emit: screened_assemblies
    path "screen_summary.tsv", emit: summary

    script:
    """
    # Reference sketch: S. aureus USA300 (NC_007793)
    REF_FASTA="${projectDir}/assets/usa300.fasta"
    REF_SKETCH="${projectDir}/assets/usa300.msh"
    
    # Download reference if not present
    if [ ! -f "$REF_FASTA" ]; then
        curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Staphylococcus_aureus/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz" | \
        gunzip -c > "$REF_FASTA"
    fi
    
    # Sketch reference if not present
    if [ ! -f "$REF_SKETCH" ]; then
        mash sketch -i "$REF_FASTA" -o "$REF_SKETCH" -k 21 -s 1000
    fi
    
    # Sketch query
    mash sketch -i "$fasta" -o "${sample_id}.msh" -k 21 -s 1000
    
    # Calculate distance
    mash dist "$REF_SKETCH" "${sample_id}.msh" > "${sample_id}.dist.tsv"
    
    # Parse and pass/fail
    DIST=$(cut -f3 "${sample_id}.dist.tsv")
    P_VALUE=$(cut -f4 "${sample_id}.dist.tsv")
    
    # Threshold: mash dist < 0.05 (~95% ANI)
    if (( $(echo "$DIST < 0.05" | bc -l) )); then
        PASS="PASS"
        ANI=$(echo "100 - ($DIST * 100)" | bc -l)
    else
        PASS="FAIL"
        ANI=$(echo "100 - ($DIST * 100)" | bc -l)
    fi
    
    # Output screen result
    echo -e "${sample_id}\\t${PASS}\\t${DIST}\\t${ANI}\\t${P_VALUE}" > "${sample_id}.screen.tsv"
    
    # Append to summary
    echo -e "${sample_id}\\t${PASS}\\t${ANI}\\t${P_VALUE}" >> screen_summary.tsv
    """
}
