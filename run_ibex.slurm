#!/bin/bash
#SBATCH --job-name=staphit
#SBATCH --output=staphit_%j.log
#SBATCH --error=staphit_%j.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=batch

module load nextflow
module load singularity

mkdir -p singularity_cache singularity_tmp
export SINGULARITY_TMPDIR=$(pwd)/singularity_tmp

# Clean previous cache to force fresh pulls
rm -rf .nextflow/cache

nextflow run main.nf \
    -profile ibex \
    --input samplesheet.csv \
    -work-dir ./work \
    -resume

