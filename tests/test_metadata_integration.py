"""Integration test: template -> fill -> validate -> normalize -> verify JSON structure."""
import subprocess
import csv
import json
import os
import tempfile
import pytest

TOOL = os.path.join(os.path.dirname(__file__), '..', 'bin', 'staphit-metadata')


def test_full_workflow():
    """End-to-end: generate template, fill it, validate, normalize."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # 1. Create samplesheet
        ss_path = os.path.join(tmpdir, 'samplesheet.csv')
        with open(ss_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['sample', 'fastq_1', 'fastq_2'])
            w.writerow(['ID00001', '/data/R1.fq.gz', '/data/R2.fq.gz'])
            w.writerow(['ID00002', '/data/R1.fq.gz', '/data/R2.fq.gz'])

        # 2. Generate templates
        result = subprocess.run(
            ['python', TOOL, 'template', '--samplesheet', ss_path,
             '--antibiogram', '--output-dir', tmpdir],
            capture_output=True, text=True
        )
        assert result.returncode == 0

        # 3. Fill in metadata
        meta_path = os.path.join(tmpdir, 'sample_metadata.csv')
        with open(meta_path) as f:
            rows = list(csv.DictReader(f))
        rows[0]['collection_date'] = '2024-06-01'
        rows[0]['geo_loc_country'] = 'Saudi Arabia'
        rows[0]['geo_loc_region'] = 'Riyadh'
        rows[0]['isolation_source'] = 'blood'
        rows[0]['infection_origin'] = 'hospital'
        rows[1]['collection_date'] = '2024-06-02'
        rows[1]['geo_loc_country'] = 'Saudi Arabia'
        rows[1]['geo_loc_region'] = 'Jeddah'
        rows[1]['isolation_source'] = 'wound'
        rows[1]['infection_origin'] = 'community'
        with open(meta_path, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=rows[0].keys())
            w.writeheader()
            w.writerows(rows)

        # 4. Fill in antibiogram
        abg_path = os.path.join(tmpdir, 'antibiogram.csv')
        with open(abg_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['sample_id', 'antibiotic', 'resistance_phenotype',
                        'measurement', 'measurement_sign', 'measurement_units',
                        'laboratory_typing_method', 'testing_standard',
                        'testing_standard_version', 'platform'])
            w.writerow(['ID00001', 'oxacillin', 'R', '4', '=', 'mg/L', 'MIC', 'CLSI', '', ''])
            w.writerow(['ID00001', 'vancomycin', 'S', '1', '=', 'mg/L', 'MIC', 'CLSI', '', ''])
            w.writerow(['ID00002', 'oxacillin', 'S', '0.5', '<=', 'mg/L', 'MIC', 'CLSI', '', ''])

        # 5. Validate
        result = subprocess.run(
            ['python', TOOL, 'validate', '--metadata', meta_path,
             '--antibiogram', abg_path, '--samplesheet', ss_path],
            capture_output=True, text=True
        )
        assert result.returncode == 0

        # 6. Normalize
        json_path = os.path.join(tmpdir, 'metadata.json')
        result = subprocess.run(
            ['python', TOOL, 'normalize', '--metadata', meta_path,
             '--antibiogram', abg_path, '-o', json_path],
            capture_output=True, text=True
        )
        assert result.returncode == 0

        # 7. Verify JSON
        with open(json_path) as f:
            data = json.load(f)

        assert len(data) == 2
        assert data[0]['run_id'] == 'ID00001'
        assert data[0]['infection_origin'] == 'hospital'
        assert len(data[0]['antibiogram']) == 2
        assert data[0]['antibiogram'][0]['antibiotic'] == 'oxacillin'
        assert data[0]['antibiogram'][0]['sir'] == 'R'
        assert data[1]['run_id'] == 'ID00002'
        assert len(data[1]['antibiogram']) == 1
