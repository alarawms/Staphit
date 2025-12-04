import argparse
import csv
import requests
import xml.etree.ElementTree as ET
import time
import os

# Base URL for NCBI E-utilities
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")

def esearch_organism(organism_name: str, max_results: int = 6) -> list:
    """
    Searches SRA for a specific organism with strict WGS/Illumina/Paired filters.
    """
    # Simplified term to ensure we get hits
    term = f'"{organism_name}"[Organism] AND "paired"[Layout] AND "illumina"[Platform]'
    
    print(f"Searching SRA for: {term}")
    params = {
        "db": "sra",
        "term": term,
        "retmode": "json",
        "retmax": 100 # Fetch more than needed to filter
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    try:
        response = requests.get(f"{EUTILS_BASE}esearch.fcgi", params=params)
        response.raise_for_status()
        data = response.json()
        return data.get("esearchresult", {}).get("idlist", [])
    except Exception as e:
        print(f"Error during search: {e}")
        return []

def efetch_details(uids: list) -> list:
    """
    Fetches XML details for a list of SRA UIDs.
    """
    if not uids:
        return []
        
    # Batch fetch
    params = {
        "db": "sra",
        "id": ",".join(uids),
        "rettype": "full",
        "retmode": "xml"
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    try:
        response = requests.get(f"{EUTILS_BASE}efetch.fcgi", params=params)
        response.raise_for_status()
        # print(f"DEBUG: XML Response Length: {len(response.text)}")
        if len(response.text) < 500:
             print(f"DEBUG: Short XML Response: {response.text}")
        
        return parse_sra_xml(response.text)
    except Exception as e:
        print(f"Error fetching details: {e}")
        return []

def parse_sra_xml(xml_string: str) -> list:
    """
    Parses SRA XML to extract run accessions and comprehensive metadata.
    """
    results = []
    try:
        # Debug: print first 500 chars of XML to see structure
        # print(f"DEBUG XML: {xml_string[:500]}") 
        
        root = ET.fromstring(xml_string)
        # print(f"DEBUG: Root tag is {root.tag}")
        
        # Check if we have any packages
        packages = root.findall(".//EXPERIMENT_PACKAGE")
        if not packages:
            print("DEBUG: No EXPERIMENT_PACKAGE found in XML.")
            return []

        for exp_pkg in packages:
            # 1. Run Info
            run_node = exp_pkg.find(".//RUN")
            if run_node is None: continue
            
            srr = run_node.get("accession")
            if not srr: continue

            # 2. Study Info
            study_node = exp_pkg.find(".//STUDY")
            study_acc = study_node.get("accession") if study_node is not None else ""
            study_title = ""
            if study_node is not None:
                desc = study_node.find("DESCRIPTOR")
                if desc is not None:
                    title_node = desc.find("STUDY_TITLE")
                    if title_node is not None:
                        study_title = title_node.text

            # 3. Sample Info & Attributes
            sample_node = exp_pkg.find(".//SAMPLE")
            sample_acc = sample_node.get("accession") if sample_node is not None else ""
            organism = "Unknown"
            if sample_node is not None:
                sample_name = sample_node.find("SAMPLE_NAME")
                if sample_name is not None:
                    # Try standard tag
                    taxon = sample_name.find("TAXON_SCIENTIFIC_NAME")
                    if taxon is not None:
                        organism = taxon.text
                    else:
                        # Try scientific_name attribute
                        sci_name = sample_name.find("SCIENTIFIC_NAME")
                        if sci_name is not None:
                            organism = sci_name.text
            
            # Dynamic Attribute Extraction
            attributes = {}
            if sample_node is not None:
                sample_attrs = sample_node.find("SAMPLE_ATTRIBUTES")
                if sample_attrs is not None:
                    for attr in sample_attrs.findall("SAMPLE_ATTRIBUTE"):
                        tag = attr.find("TAG")
                        val = attr.find("VALUE")
                        if tag is not None and val is not None and tag.text:
                            # print(f"DEBUG: Found attribute {tag.text} = {val.text}")
                            attributes[tag.text.lower().replace(" ", "_")] = val.text

            # Core metadata fields
            collection_date = attributes.get("collection_date", "not provided")
            geo_loc = attributes.get("geo_loc_name", "not provided")
            host = attributes.get("host", "not provided")
            isolation_source = attributes.get("isolation_source", "not provided")
            
            # 4. Stats (Robust extraction)
            total_spots = "0"
            try:
                # 1. Check Statistics node
                stats = run_node.find("Statistics")
                if stats is not None:
                    spots = stats.get("nspots")
                    if not spots:
                        spots_node = stats.find("Spots")
                        if spots_node is not None:
                            spots = spots_node.text
                    if spots:
                        total_spots = spots
                
                # 2. Check RUN attributes
                if total_spots == "0":
                    total_spots = run_node.get("total_spots", "0")
                    
            except:
                pass 
            
            # Filtering
            if total_spots and total_spots.isdigit() and 0 < int(total_spots) < 10000:
                # print(f"DEBUG: Skipping {srr} due to low read count: {total_spots}")
                continue 

            entry = {
                "sample": srr,
                "sra": srr,
                "fastq_1": "",
                "fastq_2": "",
                "organism": organism,
                "collection_date": collection_date,
                "geo_location": geo_loc,
                "host": host,
                "isolation_source": isolation_source,
                "study_accession": study_acc,
                "study_title": study_title,
                "read_count": total_spots,
                "library_layout": "PAIRED", 
                "platform": "ILLUMINA"
            }
            results.append(entry)
            
    except ET.ParseError as e:
        print(f"XML Parse Error: {e}")
        
    return results

def main():
    parser = argparse.ArgumentParser(description="Search SRA and generate pipeline inputs.")
    parser.add_argument("-t", "--taxid", type=str, default="Staphylococcus aureus", help="Organism name or TaxID")
    parser.add_argument("-n", "--number", type=int, default=6, help="Number of samples to fetch")
    parser.add_argument("--samplesheet", default="samplesheet.csv", help="Output samplesheet file")
    parser.add_argument("--metadata", default="metadata.csv", help="Output metadata file")
    
    args = parser.parse_args()

    # 1. Search
    uids = esearch_organism(args.taxid, max_results=args.number * 3) # Fetch extras to filter
    print(f"Found {len(uids)} potential datasets.")
    
    # 2. Fetch Details & Filter
    candidates = efetch_details(uids)
    print(f"Parsed {len(candidates)} valid candidates.")
    
    # 3. Select Top N
    selected = candidates[:args.number]
    print(f"Selected {len(selected)} samples for analysis.")

    # 4. Write Samplesheet (for Nextflow)
    ss_fields = ["sample", "sra", "fastq_1", "fastq_2"]
    with open(args.samplesheet, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=ss_fields)
        writer.writeheader()
        for s in selected:
            writer.writerow({k: s[k] for k in ss_fields})
    print(f"Written {args.samplesheet}")

    # 5. Write Metadata (Full details)
    if selected:
        md_fields = list(selected[0].keys())
        with open(args.metadata, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=md_fields)
            writer.writeheader()
            for s in selected:
                writer.writerow(s)
        print(f"Written {args.metadata}")

if __name__ == "__main__":
    main()
