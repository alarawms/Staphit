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
                ├─→ Snippy (ref-based) ─→ snippy-core ─→ IQ-TREE / SNP-dists
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
| SNP calling | Snippy | Reference-based variant calling (incremental alternative to Panaroo) |
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

# Full run with metadata and antibiogram
nextflow run main.nf -profile docker --assembler skesa \
    --metadata sample_metadata.csv \
    --antibiogram antibiogram.csv

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
| `--include_samples` | `null` | File or comma-separated list of sample IDs to include |
| `--exclude_samples` | `null` | File or comma-separated list of sample IDs to exclude |
| `--phylo_method` | `panaroo` | Phylogeny method: `panaroo` (pangenome) or `snippy` (reference-based) |
| `--reference` | `null` | Reference genome for snippy (GenBank format) |
| `--iqtree_model` | `GTR+F+I` | IQ-TREE substitution model (`MFP` for model testing) |
| `--iqtree_bb` | `1000` | Ultrafast bootstrap replicates |
| `--iqtree_fast` | `false` | Enable IQ-TREE fast mode (2-5x speedup) |
| `--iqtree_seed` | `null` | Previous `.treefile` to seed incremental tree building |

### Sample Filtering

Include or exclude samples without editing the samplesheet:

```bash
# Run only samples with metadata
cut -d, -f1 sample_metadata.csv | tail -n+2 > metadata_ids.txt
nextflow run main.nf -profile docker --include_samples metadata_ids.txt

# Exclude specific samples
nextflow run main.nf -profile docker --exclude_samples "ID00160,ID00321"

# Exclude from file (one ID per line, # comments allowed)
nextflow run main.nf -profile docker --exclude_samples bad_samples.txt
```

### Metadata and Antibiogram

#### Generating metadata from existing lab data

The `staphit-metadata convert` tool transforms various lab data formats into pipeline-ready CSVs:

```bash
# Generate sample_metadata.csv from master spreadsheet
python bin/staphit-metadata convert --from-master-csv m.csv -o sample_metadata.csv

# Generate antibiogram from Vitek 2 PDFs (full MIC data, ~560 samples)
python bin/staphit-metadata convert --from-vitek-pdf vitek_pdfs/ -o antibiogram.csv

# Generate antibiogram from wide-format Vitek TSV (SIR:MIC cells)
python bin/staphit-metadata convert --from-vitek-csv collected_metadata_res.csv -o antibiogram.csv

# Enrich metadata with clinical data from KAIMRC XLSX + extract supplementary drug rows
python bin/staphit-metadata convert \
    --from-kaimrc-xlsx "MRSA metadata and ibec KAIMRC.XLSX" \
    --metadata sample_metadata.csv \
    --existing-antibiogram antibiogram.csv \
    -o antibiogram_supplement.csv

# Merge supplementary drug rows into main antibiogram
tail -n+2 antibiogram_supplement.csv >> antibiogram.csv
```

#### Generating templates from scratch

```bash
python bin/staphit-metadata template --samplesheet samplesheet.csv --antibiogram
```

This creates empty `sample_metadata.csv` and `antibiogram.csv` templates pre-filled with sample IDs.

#### Validation

```bash
python bin/staphit-metadata validate --metadata sample_metadata.csv \
    --antibiogram antibiogram.csv --samplesheet samplesheet.csv
```

The metadata schema follows PHA4GE, MIxS v6, INSDC Pathogen.cl, and WHO GLASS standards.

### Phylogenetics

Two methods are available for building the core alignment:

#### Panaroo (default) — pangenome-based

Best for publication-quality analysis. Identifies the core genome across all samples and aligns shared genes. Must rerun from scratch when samples are added.

```bash
nextflow run main.nf -profile docker --phylo_method panaroo
```

#### Snippy — reference-based (incremental)

Best for ongoing surveillance. Maps each sample against a reference genome independently — fully cached per-sample by `-resume`. Only `snippy-core` and IQ-TREE rerun when samples are added.

```bash
# Download S. aureus NCTC 8325 reference
wget -O assets/reference.gbk \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gbff.gz"
gunzip assets/reference.gbk.gz

nextflow run main.nf -profile docker \
    --phylo_method snippy --reference assets/reference.gbk
```

#### Speeding up IQ-TREE

```bash
# Fast mode (2-5x speedup, good for surveillance)
nextflow run main.nf -profile docker --iqtree_fast true

# Seed from a previous tree (incremental, skips de novo construction)
nextflow run main.nf -profile docker \
    --iqtree_seed results/iqtree/core_gene_alignment.aln.treefile

# Full model testing (slow, publication quality)
nextflow run main.nf -profile docker --iqtree_model MFP
```

### Adding New Samples (Incremental Runs)

When new samples arrive, add them to the samplesheet and rerun with `-resume`. The pipeline caches all per-sample steps — only new samples are processed.

```bash
# 1. Add new rows to samplesheet.csv
# 2. Update metadata if needed
# 3. Rerun with resume

# Panaroo path (collective steps rerun):
nextflow run main.nf -profile docker -resume \
    --assembler skesa \
    --metadata sample_metadata.csv --antibiogram antibiogram.csv \
    --iqtree_seed results/iqtree/core_gene_alignment.aln.treefile

# Snippy path (only snippy-core + IQ-TREE rerun):
nextflow run main.nf -profile docker -resume \
    --assembler skesa --phylo_method snippy --reference assets/reference.gbk \
    --metadata sample_metadata.csv --antibiogram antibiogram.csv \
    --iqtree_seed results/iqtree/core.aln.treefile
```

**What gets cached vs rerun:**

| Step | `-resume` behavior |
|------|--------------------|
| Per-sample (Trimmomatic, SKESA, Prokka, MLST, ...) | Cached — only new samples run |
| Snippy (per-sample SNP calling) | Cached — only new samples run |
| snippy-core (merge SNPs) | Reruns (fast, minutes) |
| Panaroo (pangenome alignment) | Reruns from scratch (slow, hours) |
| IQ-TREE | Reruns, but seeded from previous tree if `--iqtree_seed` set |
| MultiQC, Summary, Visualization | Reruns (fast) |

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
├── panaroo/            # Pangenome analysis (if --phylo_method panaroo)
├── snippy/             # Per-sample SNP calls (if --phylo_method snippy)
├── snippy_core/        # Core SNP alignment (if --phylo_method snippy)
├── iqtree/             # Phylogenetic tree
├── snp_dists/          # SNP distance matrix
├── metadata/           # Validated/normalized metadata
├── aggregated/         # Per-sample summary JSONs/CSVs
├── visualization/      # Figures
└── multiqc/            # MultiQC report
```
