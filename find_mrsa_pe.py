import requests
import time
import xml.etree.ElementTree as ET

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def find_mrsa_paired():
    # Search for MRSA WGS Paired Illumina
    term = '"Staphylococcus aureus"[Organism] AND "methicillin resistant"[All Fields] AND "paired"[Layout] AND "wgs"[Strategy] AND "illumina"[Platform]'
    
    params = {
        "db": "sra",
        "term": term,
        "retmode": "json",
        "retmax": 20 
    }
    response = requests.get(f"{EUTILS_BASE}esearch.fcgi", params=params)
    data = response.json()
    id_list = data.get("esearchresult", {}).get("idlist", [])
    
    print(f"Found {len(id_list)} candidates. Checking details...")
    
    valid_srrs = []
    for uid in id_list:
        if len(valid_srrs) >= 6:
            break
            
        # Fetch summary to get SRR accession
        summary_params = {
            "db": "sra",
            "id": uid,
            "rettype": "docsum",
            "retmode": "json"
        }
        # Docsum doesn't always give SRR directly easily in JSON, let's use efetch xml for certainty or run summary
        # Actually, let's just use efetch XML like the other script
        fetch_params = {
            "db": "sra",
            "id": uid,
            "rettype": "full",
            "retmode": "xml"
        }
        r = requests.get(f"{EUTILS_BASE}efetch.fcgi", params=fetch_params)
        if r.status_code != 200: continue
        
        try:
            root = ET.fromstring(r.text)
            for run in root.findall(".//RUN"):
                srr = run.get("accession")
                # Double check layout
                layout = root.find(".//LIBRARY_LAYOUT/PAIRED")
                if layout is not None and srr:
                    valid_srrs.append(srr)
        except:
            pass
        time.sleep(0.2)

    print("Selected SRRs:", valid_srrs)
    return valid_srrs

if __name__ == "__main__":
    find_mrsa_paired()
