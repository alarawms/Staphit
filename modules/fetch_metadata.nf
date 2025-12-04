nextflow.enable.dsl=2

process FETCH_METADATA {
    container 'python:3.9-slim'
    publishDir '.', mode: 'copy', overwrite: true

    input:
    path samplesheet_csv

    output:
    path "metadata.csv", emit: metadata_csv
    path "metadata.json", emit: metadata_json

    shell:
    '''
    pip install requests

    cat << 'EOF' > fetch_metadata.py
import csv
import json
import requests
import xml.etree.ElementTree as ET
import time
import sys

# Read samples from samplesheet
runs = []
try:
    # We use the input file provided by Nextflow
    with open("!{samplesheet_csv}", "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("sra"):
                runs.append(row["sra"])
except Exception as e:
    print(f"Error reading samplesheet: {e}")
    sys.exit(1)

if not runs:
    print("No runs found in samplesheet.")
    with open("metadata.csv", "w") as f:
        print("run_id,biosample", file=f)
    with open("metadata.json", "w") as f:
        f.write("[]")
    sys.exit(0)

print(f"Fetching metadata for {len(runs)} runs...")

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# Step 1: Get RunInfo to find BioSamples
print("Step 1: getting BioSample IDs from RunInfo")
run_to_biosample = {}
biosample_to_runs = {}

for chunk in chunks(runs, 200):
    params = {
        "db": "sra",
        "id": ",".join(chunk),
        "rettype": "runinfo",
        "retmode": "text"
    }
    try:
        resp = requests.post(f"{base_url}efetch.fcgi", data=params)
        if resp.status_code != 200: continue
        
        # Use splitlines() to avoid escaping issues
        lines = resp.text.strip().splitlines()
        if len(lines) < 2: continue
        
        reader = csv.DictReader(lines)
        for row in reader:
            r = row.get("Run")
            b = row.get("BioSample")
            if r and b:
                run_to_biosample[r] = b
                if b not in biosample_to_runs: biosample_to_runs[b] = []
                biosample_to_runs[b].append(r)
    except Exception as e:
        print(f"Error in Step 1: {e}")

biosamples = list(biosample_to_runs.keys())
print(f"Found {len(biosamples)} unique BioSamples.")

# Step 2: Get Attributes from BioSample DB
# We must use ESearch to resolve Accessions (SAMEA...) to UIDs, then EFetch.
print("Step 2: getting Attributes from BioSample DB")
biosample_data = {}

for chunk in chunks(biosamples, 100):
    # Join with OR to search multiple accessions
    query = " OR ".join(chunk)
    
    # ESearch with history
    params = {
        "db": "biosample",
        "term": query,
        "usehistory": "y"
    }
    
    try:
        resp = requests.post(f"{base_url}esearch.fcgi", data=params)
        if resp.status_code != 200:
            print(f"ESearch failed: {resp.status_code}")
            continue
            
        root = ET.fromstring(resp.content)
        webenv = root.find("WebEnv")
        query_key = root.find("QueryKey")
        
        if webenv is None or query_key is None:
            print("No WebEnv found in ESearch response")
            continue
            
        webenv = webenv.text
        query_key = query_key.text
        
        # EFetch using History
        fetch_params = {
            "db": "biosample",
            "query_key": query_key,
            "WebEnv": webenv,
            "rettype": "full",
            "retmode": "xml"
        }
        
        fetch_resp = requests.post(f"{base_url}efetch.fcgi", data=fetch_params)
        if fetch_resp.status_code != 200:
            print(f"EFetch failed: {fetch_resp.status_code}")
            continue
            
        fetch_root = ET.fromstring(fetch_resp.content)
        
        for sample in fetch_root.findall("BioSample"):
            acc = sample.get("accession")
            if not acc: continue
            
            data = {}
            
            # IDs
            ids = sample.find("Ids")
            if ids is not None:
                for id_node in ids.findall("Id"):
                    if id_node.get("db") == "SRA":
                        data["sra_sample"] = id_node.text
                    if id_node.get("db") == "Sample name":
                        data["sample_name"] = id_node.text
            
            # Attributes
            attrs = sample.find("Attributes")
            if attrs is not None:
                for attr in attrs.findall("Attribute"):
                    key = attr.get("attribute_name")
                    val = attr.text
                    if key and val:
                        clean_key = attr.get("harmonized_name") or key.lower().replace(" ", "_")
                        data[clean_key] = val
            
            biosample_data[acc] = data
            
    except Exception as e:
        print(f"Error in Step 2 chunk: {e}")

# Step 3: Compile results
final_json = []
csv_rows = []
csv_headers = ["run_id", "biosample", "subject_id", "collection_date", "geo_loc_name", "host", "isolation_source", "resistance_profile"]

for run in runs:
    bs = run_to_biosample.get(run)
    if not bs: continue
    
    info = biosample_data.get(bs, {})
    
    record = {
        "run_id": run,
        "biosample": bs,
        "attributes": info
    }
    
    # Flatten for top-level JSON convenience
    if info.get("subject_id"):
        record["subject_id"] = info["subject_id"]
    for k in ["collection_date", "geo_loc_name", "host", "isolation_source"]:
        if k in info:
            record[k] = info[k]
            
    final_json.append(record)
    
    # CSV
    res_prof = info.get("antimicrobial_resistance", "") or info.get("genotype", "")
    res_prof = res_prof.replace(",", ";")
    
    row = {
        "run_id": run,
        "biosample": bs,
        "subject_id": info.get("subject_id", ""),
        "collection_date": info.get("collection_date", ""),
        "geo_loc_name": info.get("geo_loc_name", ""),
        "host": info.get("host", ""),
        "isolation_source": info.get("isolation_source", ""),
        "resistance_profile": res_prof
    }
    csv_rows.append(row)

# Write Output
with open("metadata.json", "w") as f:
    json.dump(final_json, f, indent=2)

with open("metadata.csv", "w") as f:
    writer = csv.DictWriter(f, fieldnames=csv_headers)
    writer.writeheader()
    for r in csv_rows:
        writer.writerow(r)

print("Done.")
EOF

    python fetch_metadata.py
    '''
}
