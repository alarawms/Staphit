nextflow.enable.dsl=2

process AGGREGATE {
    tag "$sample_id"
    publishDir "${params.outdir}/aggregated", mode: 'copy'
    container 'python:3.9-slim'

    input:
    tuple val(sample_id), path(trim_log), path(fastqc_files), path(quast_dir), path(mlst_tsv), path(abricate_tabs), path(metadata_json)

    output:
    path "${sample_id}_report.json"
    path "${sample_id}_summary.csv"

    script:
    """
    cat << 'EOF' > aggregate.py
import json
import csv
import os
import zipfile
import re
import sys

sample_id = "${sample_id}"
metadata_file = "${metadata_json}"
trim_log = "${trim_log}"
fastqc_files = [f for f in os.listdir('.') if f.endswith('_fastqc.zip')]
quast_dir = "${quast_dir}"
mlst_file = "${mlst_tsv}"
abricate_files = [f for f in os.listdir('.') if f.endswith('.tab') and 'abricate' in f]

data = {
    "sample_id": sample_id,
    "metadata": {},
    "qc": {},
    "assembly": {},
    "typing": {},
    "resistance": [],
    "virulence": [],
    "plasmids": []
}

# 1. Parse Metadata
try:
    if os.path.exists(metadata_file) and os.path.getsize(metadata_file) > 0:
        with open(metadata_file, 'r') as f:
            meta_list = json.load(f)
            for m in meta_list:
                if m.get("run_id") == sample_id:
                    data["metadata"] = m
                    break
except Exception as e:
    print(f"Warning: Metadata parse error: {e}")

# 2. Parse Trimmomatic
try:
    with open(trim_log, 'r') as f:
        content = f.read()
        match = re.search(r"Input Read Pairs: (\\d+).*Both Surviving: (\\d+)", content)
        if match:
            data["qc"]["raw_reads"] = int(match.group(1))
            data["qc"]["trimmed_reads"] = int(match.group(2))
            data["qc"]["survival_rate"] = (int(match.group(2)) / int(match.group(1))) * 100
except Exception as e:
    print(f"Warning: Trimmomatic parse error: {e}")

# 3. Parse FastQC
for zip_path in fastqc_files:
    try:
        with zipfile.ZipFile(zip_path, 'r') as z:
            data_file = [n for n in z.namelist() if n.endswith('fastqc_data.txt')][0]
            with z.open(data_file) as f:
                content = f.read().decode('utf-8')
                match = re.search(r"Total Sequences\\t(\\d+)", content)
                read_type = "R1" if "_1_fastqc" in zip_path else "R2"
                if match:
                    data["qc"][f"fastqc_total_{read_type}"] = int(match.group(1))
    except Exception as e:
        print(f"Warning: FastQC parse error {zip_path}: {e}")

# 4. Parse QUAST
try:
    report_path = os.path.join(quast_dir, "transposed_report.tsv")
    if os.path.exists(report_path):
        with open(report_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\\t')
            row = next(reader)
            data["assembly"]["n50"] = int(row.get("N50", 0))
            data["assembly"]["contigs"] = int(row.get("# contigs", 0))
            data["assembly"]["length"] = int(row.get("Total length", 0))
            data["assembly"]["gc"] = float(row.get("GC (%)", 0))
except Exception as e:
    print(f"Warning: QUAST parse error: {e}")

# 5. Parse MLST
try:
    with open(mlst_file, 'r') as f:
        line = f.readline()
        f.seek(0)
        delim = ',' if ',' in line else '\\t'
        reader = csv.reader(f, delimiter=delim)
        for row in reader:
            if len(row) >= 3:
                data["typing"]["mlst"] = {
                    "scheme": row[1],
                    "st": row[2],
                    "alleles": row[3:]
                }
except Exception as e:
    print(f"Warning: MLST parse error: {e}")

# 6. Parse Abricate
for tab_file in abricate_files:
    try:
        db_name = "unknown"
        if "resfinder" in tab_file: db_name = "resfinder"
        elif "vfdb" in tab_file: db_name = "vfdb"
        elif "plasmidfinder" in tab_file: db_name = "plasmidfinder"
        
        with open(tab_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\\t')
            for row in reader:
                item = {
                    "gene": row.get("GENE"),
                    "coverage": float(row.get("%COVERAGE", 0)),
                    "identity": float(row.get("%IDENTITY", 0)),
                    "accession": row.get("ACCESSION"),
                    "product": row.get("PRODUCT"),
                    "database": db_name
                }
                
                if db_name == "resfinder":
                    data["resistance"].append(item)
                elif db_name == "vfdb":
                    data["virulence"].append(item)
                elif db_name == "plasmidfinder":
                    data["plasmids"].append(item)
                    
    except Exception as e:
        print(f"Warning: Abricate parse error {tab_file}: {e}")

# 7. Output
with open(f"{sample_id}_report.json", 'w') as f:
    json.dump(data, f, indent=2)

with open(f"{sample_id}_summary.csv", 'w') as f:
    writer = csv.writer(f)
    header = ["sample_id", "st", "n50", "contigs", "length", "resistance_genes", "virulence_genes"]
    writer.writerow(header)
    
    st = data["typing"].get("mlst", {}).get("st", "ND")
    n50 = data["assembly"].get("n50", 0)
    contigs = data["assembly"].get("contigs", 0)
    length = data["assembly"].get("length", 0)
    res_genes = ";".join([r["gene"] for r in data["resistance"]])
    vir_genes = ";".join([v["gene"] for v in data["virulence"]])
    
    writer.writerow([sample_id, st, n50, contigs, length, res_genes, vir_genes])

print(f"Aggregation complete for {sample_id}")
EOF

    python aggregate.py
    """)
}