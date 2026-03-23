nextflow.enable.dsl=2

process GENERATE_SUMMARY {
    tag "summary"
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.html"

    input:
    path screen_summary
    path qc_summary
    path quast_dir
    tuple val(sample_id), path(mlst_tsv)
    tuple val(sample_id), path(spatyper_tsv)
    tuple val(sample_id), path(sccmec_tsv)
    tuple val(sample_id), path(agr_tsv)
    tuple val(sample_id), path(abricate_tsv)
    tuple val(sample_id), path(amrfinder_tsv)

    output:
    path "summary.html", emit: html
    path "summary_data.json", emit: json

    script:
    """
    # Collect all data into JSON
    python3 << 'PYEOF'
import json
import os
from pathlib import Path

# Base directories
quast_dir = "quast" if os.path.isdir("${quast_dir}") else None

# Summary data
summary = {
    "screen_species": None,
    "qc_filter": None,
    "samples": []
}

# Parse species screen summary
if os.path.exists("${screen_summary}") and os.path.getsize("${screen_summary}") > 0:
    with open("${screen_summary}") as f:
        qc_lines = [l for l in f if l.strip() and not l.startswith("sample")]
        summary["screen_species"] = {
            "total": len(qc_lines),
            "pass": len([l for l in qc_lines if "PASS" in l]),
            "fail": len([l for l in qc_lines if "FAIL" in l])
        }

# Parse QC filter summary
if os.path.exists("${qc_summary}") and os.path.getsize("${qc_summary}") > 0:
    with open("${qc_summary}") as f:
        qc_lines = [l for l in f if l.strip() and not l.startswith("sample")]
        summary["qc_filter"] = {
            "total": len(qc_lines),
            "pass": len([l for l in qc_lines if "PASS" in l]),
            "fail": len([l for l in qc_lines if "FAIL" in l])
        }

# Collect sample data
samples = {}
all_sample_ids = set()

# Collect from all typing results
tsv_files = [
    ("${mlst_tsv}", "mlst"),
    ("${spatyper_tsv}", "spatyper"),
    ("${sccmec_tsv}", "sccmec"),
    ("${agr_tsv}", "agr"),
    ("${abricate_tsv}", "abricate"),
    ("${amrfinder_tsv}", "amrfinder")
]

for tsv_file, key in tsv_files:
    if os.path.exists(tsv_file) and os.path.getsize(tsv_file) > 0:
        with open(tsv_file) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    parts = line.strip().split("\\t")
                    sample_id = parts[0]
                    all_sample_ids.add(sample_id)
                    if sample_id not in samples:
                        samples[sample_id] = {}
                    samples[sample_id][key] = "\\t".join(parts[1:])

# Parse QUAST metrics
if quast_dir:
    for sample_id in all_sample_ids:
        sample_quast = os.path.join(quast_dir, f"{sample_id}_quast_report/report.txt")
        if os.path.exists(sample_quast):
            with open(sample_quast) as f:
                for line in f:
                    if "Total length" in line:
                        samples[sample_id]["total_length"] = line.split("\\t")[1]
                    elif "# contigs" in line:
                        samples[sample_id]["contigs"] = line.split("\\t")[1]
                    elif "N50" in line:
                        samples[sample_id]["n50"] = line.split("\\t")[1]

summary["samples"] = sorted(
    [{"id": sid, **data} for sid, data in samples.items()],
    key=lambda x: x["id"]
)

# Write JSON
with open("summary_data.json", "w") as f:
    json.dump(summary, f, indent=2)

print(f"Summary generated for {len(summary['samples'])} samples")
PYEOF

# Render HTML with Jinja2
python3 << 'HTMLEOF'
from jinja2 import Template
import json

# Load JSON data
with open("summary_data.json") as f:
    data = json.load(f)

# Load and render template
with open("${projectDir}/templates/summary.html.j2") as f:
    template = Template(f.read())

html = template.render(
    samples=data["samples"],
    screen_species=data.get("screen_species"),
    qc_filter=data.get("qc_filter")
)

with open("summary.html", "w") as f:
    f.write(html)

print("HTML summary generated: summary.html")
HTMLEOF
    """
}
