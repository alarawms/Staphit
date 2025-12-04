nextflow.enable.dsl=2

process AGGREGATOR {
    tag "$sample_id"
    publishDir "${params.outdir}/aggregated", mode: 'copy'
    container 'python:3.9-slim'

    input:
    tuple val(sample_id), path(fastp_json), path(quast_dir), path(mlst_tsv), path(abricate_tabs), path(amrfinder_report), path(mash_sketch), path(spa_report), path(sccmec_report), path(agr_report), path(kma_res), path(metadata_json)

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
fastp_file = "${fastp_json}"
quast_dir = "${quast_dir}"
mlst_file = "${mlst_tsv}"
abricate_files = [f for f in os.listdir('.') if f.endswith('.tab') and 'abricate' in f]
amrfinder_file = "${amrfinder_report}"
mash_file = "${mash_sketch}"
spa_file = "${spa_report}"
sccmec_file = "${sccmec_report}"
agr_file = "${agr_report}"
kma_file = "${kma_res}"

data = {
    "sample_id": sample_id,
    "metadata": {},
    "qc": {},
    "assembly": {},
    "typing": {
        "mlst": {},
        "spa": {},
        "sccmec": {},
        "agr": {}
    },
    "resistance": {
        "abricate": [],
        "amrfinder": [],
        "kma": []
    },
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

# 2. Parse Fastp (QC)
try:
    with open(fastp_file, 'r') as f:
        fastp_data = json.load(f)
        data["qc"]["raw_reads"] = fastp_data.get("summary", {}).get("before_filtering", {}).get("total_reads", 0)
        data["qc"]["trimmed_reads"] = fastp_data.get("summary", {}).get("after_filtering", {}).get("total_reads", 0)
        if data["qc"]["raw_reads"] > 0:
            data["qc"]["survival_rate"] = (data["qc"]["trimmed_reads"] / data["qc"]["raw_reads"]) * 100
        else:
            data["qc"]["survival_rate"] = 0
        data["qc"]["q30_rate"] = fastp_data.get("summary", {}).get("after_filtering", {}).get("q30_rate", 0)
except Exception as e:
    print(f"Warning: Fastp parse error: {e}")

# 3. Parse QUAST
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

# 4. Parse MLST
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

# 5. Parse SpaTyper
try:
    with open(spa_file, 'r') as f:
        reader = csv.reader(f, delimiter='\\t')
        header = next(reader, None)
        row = None
        if header:
            if "Repeats" in header or "Type" in header:
                row = next(reader, None)
            else:
                row = header
        if row and len(row) >= 3:
             data["typing"]["spa"] = {
                 "type": row[2],
                 "repeats": row[1]
             }
except Exception as e:
    print(f"Warning: SpaTyper parse error: {e}")

# 6. Parse SCCmec
try:
    with open(sccmec_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\\t')
        row = next(reader, None)
        if row:
             sccmec_type = row.get("SCCmec_Type") or row.get("type") or "ND"
             data["typing"]["sccmec"] = {
                 "type": sccmec_type,
                 "full_row": row
             }
except Exception as e:
    print(f"Warning: SCCmec parse error: {e}")

# 7. Parse AgrVATE (Now Staph Agr Typer JSON)
try:
    if os.path.exists(agr_file) and os.path.getsize(agr_file) > 0:
        with open(agr_file, 'r') as f:
            agr_data = json.load(f)
            data["typing"]["agr"] = {
                "group": agr_data.get("agr_group", "ND"),
                "confidence": agr_data.get("confidence", 0.0)
            }
except Exception as e:
    print(f"Warning: AgrVATE parse error: {e}")

# 8. Parse Abricate
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
                    data["resistance"]["abricate"].append(item)
                elif db_name == "vfdb":
                    data["virulence"].append(item)
                elif db_name == "plasmidfinder":
                    data["plasmids"].append(item)
                    
    except Exception as e:
        print(f"Warning: Abricate parse error {tab_file}: {e}")

# 9. Parse AMRFinderPlus
try:
    with open(amrfinder_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\\t')
        for row in reader:
            item = {
                "gene": row.get("Gene symbol"),
                "coverage": float(row.get("% Coverage", 0)),
                "identity": float(row.get("% Identity", 0)),
                "accession": row.get("Accession of closest sequence"),
                "product": row.get("Sequence name"),
                "element_type": row.get("Element type"),
                "element_subtype": row.get("Element subtype"),
                "class": row.get("Class"),
                "subclass": row.get("Subclass")
            }
            data["resistance"]["amrfinder"].append(item)
except Exception as e:
    print(f"Warning: AMRFinderPlus parse error: {e}")

# 10. Parse KMA (ResFinder)
try:
    if os.path.exists(kma_file):
        with open(kma_file, 'r') as f:
            # Skip header lines usually starting with #
            # KMA output format: Template, Score, Expected, Template_length, Template_Identity, ...
            reader = csv.reader(f, delimiter='\\t')
            header = next(reader, None)
            while header and header[0].startswith('#'):
                 header = next(reader, None)
            
            # If header is valid, process rows
            if header:
                for row in reader:
                    if len(row) > 0:
                        item = {
                            "gene": row[0], # Template name
                            "score": row[1],
                            "identity": row[4]
                        }
                        data["resistance"]["kma"].append(item)
except Exception as e:
    print(f"Warning: KMA parse error: {e}")


# 11. Output
with open(f"{sample_id}_report.json", 'w') as f:
    json.dump(data, f, indent=2)

# Generate Bactopia-style summary CSV
with open(f"{sample_id}_summary.csv", 'w') as f:
    writer = csv.writer(f, delimiter='\\t')
    
    # Define columns similar to Bactopia
    headers = [
        "sample_id", 
        "total_reads", "trimmed_reads", "survival_rate", "q30_rate", # QC
        "assembly_length", "contigs", "n50", "gc_percent", # Assembly
        "mlst_scheme", "mlst_st", # Typing
        "spa_type", "sccmec_type", "agr_group", # MRSA Typing
        "amrfinder_genes", "abricate_genes", "kma_genes", # AMR
        "virulence_genes", "plasmids" # Other
    ]
    writer.writerow(headers)
    
    # Extract values
    total_reads = data["qc"].get("raw_reads", 0)
    trimmed_reads = data["qc"].get("trimmed_reads", 0)
    survival_rate = f"{data['qc'].get('survival_rate', 0):.2f}"
    q30_rate = f"{data['qc'].get('q30_rate', 0):.2f}"
    
    asm_len = data["assembly"].get("length", 0)
    contigs = data["assembly"].get("contigs", 0)
    n50 = data["assembly"].get("n50", 0)
    gc = f"{data['assembly'].get('gc', 0):.2f}"
    
    mlst_scheme = data["typing"].get("mlst", {}).get("scheme", "-")
    mlst_st = data["typing"].get("mlst", {}).get("st", "-")
    
    spa_type = data["typing"].get("spa", {}).get("type", "-")
    sccmec_type = data["typing"].get("sccmec", {}).get("type", "-")
    agr_group = data["typing"].get("agr", {}).get("group", "-")
    
    amr_genes = ";".join(sorted(list(set([r["gene"] for r in data["resistance"]["amrfinder"] if r["gene"]]))))
    abricate_genes = ";".join(sorted(list(set([r["gene"] for r in data["resistance"]["abricate"] if r["gene"]]))))
    kma_genes = ";".join(sorted(list(set([r["gene"] for r in data["resistance"]["kma"] if r["gene"]]))))
    vir_genes = ";".join(sorted(list(set([v["gene"] for v in data["virulence"] if v["gene"]]))))
    plasmids = ";".join(sorted(list(set([p["gene"] for p in data["plasmids"] if p["gene"]]))))
    
    writer.writerow([
        sample_id,
        total_reads, trimmed_reads, survival_rate, q30_rate,
        asm_len, contigs, n50, gc,
        mlst_scheme, mlst_st,
        spa_type, sccmec_type, agr_group,
        amr_genes, abricate_genes, kma_genes,
        vir_genes, plasmids
    ])

print(f"Aggregation complete for {sample_id}")
EOF

    python aggregate.py
    """
}
