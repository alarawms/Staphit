import argparse
import csv
import requests
import xml.etree.ElementTree as ET
import os
import re
import time

# Base URL for NCBI E-utilities
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_API_KEY = os.environ.get("NCBI_API_KEY") # Optional: set NCBI_API_KEY environment variable

def esearch_sra(term: str) -> list:
    """
    Searches the SRA database for a given term (SRR, SRP, SRX) and returns a list of UIDs.
    If an SRP is given, it will return UIDs for all associated SRR runs.
    """
    params = {
        "db": "sra",
        "term": term,
        "retmode": "json",
        "retmax": 10000 # Max UIDs
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    try:
        response = requests.get(f"{EUTILS_BASE}esearch.fcgi", params=params)
        response.raise_for_status()
        data = response.json()
        return data.get("esearchresult", {}).get("idlist", [])
    except requests.exceptions.RequestException as e:
        print(f"Error during ESearch for term '{term}': {e}")
        return []

def efetch_sra_xml(uid: str) -> str | None:
    """
    Fetches detailed SRA metadata in XML format for a given UID.
    """
    params = {
        "db": "sra",
        "id": uid,
        "rettype": "full",
        "retmode": "xml"
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    try:
        response = requests.get(f"{EUTILS_BASE}efetch.fcgi", params=params)
        response.raise_for_status()
        return response.text
    except requests.exceptions.RequestException as e:
        print(f"Error during EFetch for UID '{uid}': {e}")
        return None

def parse_sra_xml_to_samples(xml_string: str) -> list:
    """
    Parses SRA XML string (from EFetch) and extracts metadata for each run.
    Returns a list of dictionaries, one for each SRR found.
    """
    samples_data = []
    try:
        root = ET.fromstring(xml_string)

        for experiment_package in root.findall(".//EXPERIMENT_PACKAGE"):
            study_accession = experiment_package.find(".//STUDY_REF").get("accession") if experiment_package.find(".//STUDY_REF") is not None else ""
            study_title = experiment_package.find(".//STUDY/DESCRIPTOR/STUDY_TITLE").text if experiment_package.find(".//STUDY/DESCRIPTOR/STUDY_TITLE") is not None else ""

            sample_accession = experiment_package.find(".//SAMPLE_REF").get("accession") if experiment_package.find(".//SAMPLE_REF") is not None else ""
            sample_title = experiment_package.find(".//SAMPLE/TITLE").text if experiment_package.find(".//SAMPLE/TITLE") is not None else ""
            organism = experiment_package.find(".//SAMPLE/SAMPLE_NAME/TAXON_SCIENTIFIC_NAME").text if experiment_package.find(".//SAMPLE/SAMPLE_NAME/TAXON_SCIENTIFIC_NAME") is not None else ""
            
            # Extract attributes like collection_date, geo_loc_name, host, isolation_source
            collection_date = ""
            geo_location = ""
            host = ""
            isolation_source = ""
            
            for attr in experiment_package.findall(".//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
                tag = attr.find("TAG")
                value = attr.find("VALUE")
                if tag is not None and value is not None:
                    if tag.text == "collection_date":
                        collection_date = value.text
                    elif tag.text == "geo_loc_name":
                        geo_location = value.text
                    elif tag.text == "host":
                        host = value.text
                    elif tag.text == "isolation_source":
                        isolation_source = value.text
            
            for run_node in experiment_package.findall(".//RUN"):
                srr_accession = run_node.get("accession")
                library_layout_node = run_node.find(".//LIBRARY_LAYOUT")
                is_paired = library_layout_node.find("PAIRED") is not None if library_layout_node is not None else False
                total_spots = run_node.find(".//Statistics/Spots").text if run_node.find(".//Statistics/Spots") is not None else ""

                samples_data.append({
                    "sample": srr_accession, # Use SRR as sample name for SRA downloads
                    "sra": srr_accession,
                    "fastq_1": "",
                    "fastq_2": "",
                    "organism": organism,
                    "collection_date": collection_date,
                    "geo_location": geo_location,
                    "host": host,
                    "isolation_source": isolation_source,
                    "study_accession": study_accession,
                    "study_title": study_title,
                    "read_count_raw": total_spots,
                    "library_layout": "PAIRED" if is_paired else "SINGLE"
                })
        
    except ET.ParseError as e:
        print(f"Error parsing SRA XML: {e}")
    return samples_data

def process_sra_accession(accession: str) -> list:
    """
    Fetches and processes metadata for a single SRA accession (SRR, SRP, SRX).
    Returns a list of sample dictionaries.
    """
    print(f"Processing SRA accession: {accession}")
    uids = esearch_sra(accession)
    if not uids:
        print(f"Warning: No SRA UIDs found for '{accession}'. Skipping.")
        return []
    
    all_sra_metadata = []
    for uid in uids:
        xml_data = efetch_sra_xml(uid)
        if xml_data:
            all_sra_metadata.extend(parse_sra_xml_to_samples(xml_data))
        else:
            print(f"Warning: Could not fetch XML for UID '{uid}'. Skipping.")
        time.sleep(0.1) # Be kind to NCBI servers
    
    # Filter to unique SRR accessions
    unique_srr_samples = {s["sample"]: s for s in all_sra_metadata if s["sra"].startswith("SRR")}.values()
    return list(unique_srr_samples)

def process_local_files(r1_path: str, r2_path: str = None) -> dict:
    """
    Processes local FastQ files and returns a sample dictionary.
    """
    if not os.path.exists(r1_path):
        raise FileNotFoundError(f"Local FastQ R1 not found: {r1_path}")
    if r2_path and not os.path.exists(r2_path):
        raise FileNotFoundError(f"Local FastQ R2 not found: {r2_path}")

    # Generate a sample name from R1 filename
    sample_name = os.path.basename(r1_path)
    # Remove common FastQ extensions and paired-end indicators
    sample_name = re.sub(r'(_R1|_1)?(_001)?(\.fastq|\.fq)?(\.gz)?$', '', sample_name, flags=re.IGNORECASE)
    
    return {
        "sample": sample_name,
        "sra": "",
        "fastq_1": os.path.abspath(r1_path),
        "fastq_2": os.path.abspath(r2_path) if r2_path else "",
        "organism": "unknown",
        "collection_date": "",
        "geo_location": "",
        "study_accession": "",
        "study_title": "",
        "read_count_raw": "",
        "library_layout": "PAIRED" if r2_path else "SINGLE"
    }

def find_fastq_pairs(directory: str) -> list:
    """
    Scans a directory for FastQ files and pairs them based on common suffixes.
    Returns a list of (r1, r2) tuples.
    """
    fastq_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if re.search(r'\.(fastq|fq)(\.gz)?$', file, re.IGNORECASE):
                fastq_files.append(os.path.join(root, file))
    
    pairs = {}
    # Common paired-end suffixes
    # (suffix_pattern, r1_marker, r2_marker)
    # We try to match R1 first, then look for R2
    patterns = [
        (r'(_R1)(_001)?\.(fastq|fq)(\.gz)?$', '_R1', '_R2'),
        (r'(_1)(_001)?\.(fastq|fq)(\.gz)?$', '_1', '_2')
    ]

    processed_files = set()

    for f in sorted(fastq_files):
        if f in processed_files:
            continue
            
        matched = False
        for pattern, r1_marker, r2_marker in patterns:
            match = re.search(pattern, f, re.IGNORECASE)
            if match:
                # This is a potential R1 file
                r1_path = f
                # Construct expected R2 path
                r2_path = f.replace(r1_marker, r2_marker)
                
                if os.path.exists(r2_path) and r2_path in fastq_files:
                    pairs[r1_path] = r2_path
                    processed_files.add(r1_path)
                    processed_files.add(r2_path)
                    matched = True
                    break
        
        if not matched:
            # If it didn't match R1 pattern, check if it's an R2 file that we haven't processed yet
            # (though sorted order usually handles R1 first).
            # Or it might be a single-end file (which we currently skip or handle as single?)
            # For now, we only return pairs as requested.
            pass

    return list(pairs.items())

def main():
    parser = argparse.ArgumentParser(
        description="Generate a samplesheet.csv for the MRSA Nextflow pipeline.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-o", "--output", default="samplesheet_prepared.csv",
                        help="Output samplesheet CSV file name (default: samplesheet_prepared.csv)")
    parser.add_argument("-s", "--sra", nargs='*', default=[],
                        help="List of SRA accessions (SRR, SRP, SRX) to include. "
                             "e.g., -s SRR12142664 SRR12142665 SRP036483")
    parser.add_argument("-f", "--fastq", nargs='*', default=[],
                        help="List of local FastQ file pairs (R1 R2). "
                             "Provide R1 and R2 paths separated by a space. "
                             "e.g., -f data/sampleA_R1.fastq.gz data/sampleA_R2.fastq.gz")
    parser.add_argument("-d", "--input-dir", 
                        help="Directory to scan for paired FastQ files. "
                             "Automatically pairs files with _R1/_R2 or _1/_2 suffixes.")
    
    args = parser.parse_args()

    all_samples = []

    # Process SRA accessions
    for sra_accession in args.sra:
        sra_records = process_sra_accession(sra_accession)
        all_samples.extend(sra_records)

    # Process local FastQ files from manual list
    if args.fastq:
        if len(args.fastq) % 2 != 0:
            print("Error: Local FastQ files must be provided in R1 R2 pairs.")
            exit(1)
        
        for i in range(0, len(args.fastq), 2):
            r1 = args.fastq[i]
            r2 = args.fastq[i+1]
            try:
                local_record = process_local_files(r1, r2)
                all_samples.append(local_record)
            except FileNotFoundError as e:
                print(f"Error processing local files: {e}. Skipping.")
            except Exception as e:
                print(f"An unexpected error occurred for local files {r1}, {r2}: {e}. Skipping.")

    # Process local FastQ files from directory
    if args.input_dir:
        if not os.path.isdir(args.input_dir):
            print(f"Error: Input directory '{args.input_dir}' does not exist.")
        else:
            print(f"Scanning directory '{args.input_dir}' for FastQ pairs...")
            pairs = find_fastq_pairs(args.input_dir)
            if not pairs:
                print("No paired FastQ files found in directory.")
            
            for r1, r2 in pairs:
                try:
                    local_record = process_local_files(r1, r2)
                    all_samples.append(local_record)
                    print(f"Found pair: {os.path.basename(r1)} / {os.path.basename(r2)}")
                except Exception as e:
                    print(f"Error processing pair {r1}, {r2}: {e}")

    if not all_samples:
        print("No samples processed. Samplesheet will not be created.")
        return

    # Define CSV fieldnames
    fieldnames = [
        "sample", "sra", "fastq_1", "fastq_2", "organism",
        "collection_date", "geo_location", "host", "isolation_source",
        "study_accession", "study_title", "read_count_raw", "library_layout"
    ]

    # Write samplesheet
    with open(args.output, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for sample_data in all_samples:
            writer.writerow(sample_data)
    
    print(f"\nSuccessfully created samplesheet: {args.output} with {len(all_samples)} entries.")

if __name__ == "__main__":
    main()
