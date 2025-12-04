# Staphit - MRSA Pathogen Tracker

A Nextflow pipeline for MRSA genomic surveillance.

## Introduction

This pipeline performs a series of analyses on MRSA (Methicillin-resistant *Staphylococcus aureus*) genomic data. It takes raw sequencing reads as input and performs the following steps:

1.  **Data Fetching:** Downloads raw sequencing data from the Sequence Read Archive (SRA).
2.  **Quality Control:** Trims low-quality bases and adapters using Trimmomatic and assesses read quality with FastQC.
3.  **Assembly:** Assembles the trimmed reads into a draft genome using SPAdes.
4.  **Annotation:** Annotates the assembled genome using Prokka.
5.  **AMR Gene Detection:** Identifies antimicrobial resistance genes using ABRicate.
6.  **MLST Typing:** Determines the Multi-Locus Sequence Type (MLST) of the isolate.
7.  **Reporting:** Aggregates all results into a single report using MultiQC.

## Dependencies

This pipeline requires [Nextflow](https://www.nextflow.io/) and [Docker](https://www.docker.com/) to be installed.

## Usage

### 1. Create a sample sheet

Create a CSV file named `samplesheet.csv` with two columns: `sample` and `sra`.

**Example `samplesheet.csv`:**
```csv
sample,sra
MRSA_TEST_1,SRR935090
MRSA_TEST_2,SRR935091
```

### 2. Run the pipeline

Run the pipeline using the following command:

```bash
nextflow run main.nf -profile docker
```

You can specify a different sample sheet using the `--input` parameter:
```bash
nextflow run main.nf -profile docker --input /path/to/your/samplesheet.csv
```

## Output

The pipeline will create a `results` directory containing the output of each step. The final report will be located in `results/multiqc/multiqc_report.html`.
