# Staphit - MRSA Genomic Surveillance Pipeline

A Nextflow (DSL2) pipeline for comprehensive MRSA (*Staphylococcus aureus*) genomic surveillance — from raw reads to resistance profiles, typing, phylogenetics, and aggregated reports.

## Pipeline Overview

```
Reads ─→ Trimmomatic ─→ FastQC
                │
                ├─→ SPAdes / SKESA ─→ Assembly QC (QUAST)
                │         │
                │         ├─→ Prokka ─→ Panaroo ─→ IQ-TREE / SNP-dists
                │         ├─→ ABRicate
                │         ├─→ AMRFinderPlus
                │         ├─→ MLST
                │         ├─→ Mash
                │         ├─→ spaTyper
                │         ├─→ SCCmec Typer
                │         └─→ agr Typing
                │
                └─→ KMA (ResFinder DB)
                          │
                          └─→ Aggregator ─→ Summary ─→ Visualization
                                                     ─→ MultiQC
```

### Analysis Steps

| Step | Tool | Description |
|------|------|-------------|
| Data fetching | SRA Tools | Download reads from NCBI SRA (or use local FASTQs) |
| Read QC | Trimmomatic + FastQC | Adapter trimming and quality assessment |
| Assembly | SPAdes or SKESA | *De novo* genome assembly (selectable via `--assembler`) |
| Assembly QC | QUAST | Assembly quality metrics |
| Species ID | Mash | Genomic distance-based species confirmation |
| Annotation | Prokka | Gene prediction and functional annotation |
| AMR detection | ABRicate + AMRFinderPlus + KMA | Resistance gene identification (assembly- and read-based) |
| MLST | mlst | Multi-locus sequence typing |
| spa typing | spaTyper | *spa* gene repeat typing |
| SCCmec typing | SCCmec Typer | Staphylococcal cassette chromosome *mec* classification |
| agr typing | agr Typer | Accessory gene regulator group assignment |
| Pangenome | Panaroo | Core/accessory genome analysis |
| Phylogenetics | IQ-TREE + SNP-dists | Maximum-likelihood tree and pairwise SNP distances |
| Reporting | MultiQC + Aggregator | Per-sample summaries and combined report |

## Requirements

- [Nextflow](https://www.nextflow.io/) (>= 22.10)
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/)

## Usage

### Input

The pipeline accepts a CSV samplesheet with either SRA accessions or local FASTQ paths:

```csv
sample,sra
MRSA_001,SRR935090
MRSA_002,SRR935091
```

```csv
sample,fastq_1,fastq_2
MRSA_001,/path/to/reads_R1.fastq.gz,/path/to/reads_R2.fastq.gz
MRSA_002,/path/to/reads_R1.fastq.gz,/path/to/reads_R2.fastq.gz
```

### Running the Pipeline

```bash
# Basic run with SPAdes (default)
nextflow run main.nf -profile docker

# Use SKESA assembler
nextflow run main.nf -profile docker --assembler skesa

# Custom samplesheet
nextflow run main.nf -profile docker --input samplesheet.csv

# Auto-search SRA for a species
nextflow run main.nf -profile docker --species "Staphylococcus aureus" --search_limit 200

# Resume a previous run
nextflow run main.nf -profile docker -resume
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | `samplesheet.csv` | Path to input samplesheet |
| `--outdir` | `./results` | Output directory |
| `--assembler` | `spades` | Assembler to use: `spades` or `skesa` |
| `--species` | `null` | Species name for SRA search (skips samplesheet) |
| `--search_limit` | `100` | Max SRA records to return |
| `--max_downloads` | `50` | Max SRA samples to download |
| `--reads_limit` | `null` | Limit number of reads per sample |
| `--metadata` | `null` | Path to sample metadata CSV |
| `--antibiogram` | `null` | Path to antibiogram CSV (long format) |

### Metadata Input

Generate a metadata template pre-filled with your sample IDs:

```bash
python bin/staphit-metadata template --samplesheet samplesheet.csv --antibiogram
```

This creates `sample_metadata.csv` and `antibiogram.csv`. Fill them in and pass to the pipeline:

```bash
nextflow run main.nf -profile docker --input samplesheet.csv \
    --metadata sample_metadata.csv \
    --antibiogram antibiogram.csv
```

Validate before running (optional):

```bash
python bin/staphit-metadata validate --metadata sample_metadata.csv \
    --antibiogram antibiogram.csv --samplesheet samplesheet.csv
```

The metadata schema follows PHA4GE, MIxS v6, INSDC Pathogen.cl, and WHO GLASS standards. See `docs/plans/2026-03-24-metadata-feature-design.md` for full field definitions.

### Profiles

| Profile | Description |
|---------|-------------|
| `docker` | Run with Docker containers (local workstation) |
| `ibex` | Run on KAUST Ibex HPC cluster with SLURM + Singularity |

## Output

Results are organized under `./results/`:

```
results/
├── trimmomatic/        # Trimmed reads
├── fastqc/             # Read quality reports
├── spades/ or skesa/   # Assembled genomes
├── quast/              # Assembly QC
├── prokka/             # Genome annotations
├── abricate/           # ABRicate AMR results
├── amrfinderplus/      # AMRFinderPlus results
├── kma/                # KMA read-based AMR
├── mlst/               # MLST types
├── mash/               # Mash distance sketches
├── spatyper/           # spa types
├── sccmec/             # SCCmec types
├── agr_typing/         # agr groups
├── panaroo/            # Pangenome analysis
├── iqtree/             # Phylogenetic tree
├── snp_dists/          # SNP distance matrix
├── metadata/           # Validated/normalized metadata
├── aggregated/         # Per-sample summary JSONs/CSVs
├── visualization/      # Figures
└── multiqc/            # MultiQC report
```
