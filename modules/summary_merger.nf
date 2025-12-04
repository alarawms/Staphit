nextflow.enable.dsl=2

process SUMMARY_MERGER {
    publishDir "${params.outdir}/aggregated", mode: 'copy'
    container 'python:3.9-slim'

    input:
    path summary_files

    output:
    path "final_summary.tsv"

    script:
    """
    cat << 'EOF' > merge_summaries.py
import csv
import sys
import os

output_file = "final_summary.tsv"
input_files = "${summary_files}".split()

if not input_files:
    print("No summary files to merge.")
    sys.exit(0)

# Read headers from the first file
with open(input_files[0], 'r') as f:
    reader = csv.reader(f, delimiter='\\t')
    header = next(reader)

with open(output_file, 'w', newline='') as fout:
    writer = csv.writer(fout, delimiter='\\t')
    writer.writerow(header)

    for fname in input_files:
        with open(fname, 'r') as fin:
            reader = csv.reader(fin, delimiter='\\t')
            try:
                file_header = next(reader)
                # Verify header matches (optional, but good practice)
                if file_header != header:
                    print(f"Warning: Header mismatch in {fname}. Skipping.")
                    continue
                
                for row in reader:
                    writer.writerow(row)
            except StopIteration:
                print(f"Warning: Empty file {fname}")
                continue

print(f"Merged {len(input_files)} summary files into {output_file}")
EOF

    python merge_summaries.py
    """
}
