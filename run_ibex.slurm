#!/bin/bash
#SBATCH --job-name=staphit
#SBATCH --output=staphit_%j.log
#SBATCH --error=staphit_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=batch

# ============================================
# Staphit Pipeline - IBEX Submission Script
# ============================================

module load nextflow
module load singularity

# Set environment
export NXF_SINGULARITY_CACHEDIR=$(pwd)/singularity_cache
export SINGULARITY_CACHEDIR=$(pwd)/singularity_cache
export SINGULARITY_TMPDIR=$(pwd)/singularity_tmp

# Create directories
mkdir -p $SINGULARITY_CACHEDIR $SINGULARITY_TMPDIR

# Check if containers exist
CONTAINER_COUNT=$(ls -1 $SINGULARITY_CACHEDIR/*.sif 2>/dev/null | wc -l)
echo "Found $CONTAINER_COUNT cached containers"

if [ "$CONTAINER_COUNT" -lt 5 ]; then
    echo "WARNING: Few containers cached. Consider running pull_containers.slurm first"
    echo "Continuing anyway - Nextflow will try to pull missing containers..."
fi

# Run the pipeline
nextflow run main.nf \
    -profile ibex \
    --input samplesheet.csv \
    -work-dir ./work \
    -resume

echo "Pipeline completed at $(date)"

