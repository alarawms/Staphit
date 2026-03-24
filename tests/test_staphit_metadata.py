"""Tests for staphit-metadata CLI tool."""
import subprocess
import csv
import os
import tempfile
import pytest

TOOL = os.path.join(os.path.dirname(__file__), '..', 'bin', 'staphit-metadata')

@pytest.fixture
def tmpdir():
    with tempfile.TemporaryDirectory() as d:
        yield d

@pytest.fixture
def samplesheet(tmpdir):
    path = os.path.join(tmpdir, 'samplesheet.csv')
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['sample', 'fastq_1', 'fastq_2'])
        writer.writerow(['ID00001', '/data/ID00001_R1.fastq.gz', '/data/ID00001_R2.fastq.gz'])
        writer.writerow(['ID00002', '/data/ID00002_R1.fastq.gz', '/data/ID00002_R2.fastq.gz'])
        writer.writerow(['ID00003', '/data/ID00003_R1.fastq.gz', '/data/ID00003_R2.fastq.gz'])
    return path


class TestTemplate:
    def test_blank_template(self, tmpdir):
        """Generate blank metadata template with header only."""
        out = os.path.join(tmpdir, 'sample_metadata.csv')
        result = subprocess.run(
            ['python', TOOL, 'template', '--output-dir', tmpdir],
            capture_output=True, text=True
        )
        assert result.returncode == 0
        assert os.path.exists(out)
        with open(out) as f:
            reader = csv.reader(f)
            header = next(reader)
            assert header[0] == 'sample_id'
            assert 'organism' in header
            assert 'collection_date' in header
            rows = list(reader)
            assert len(rows) == 0  # blank, no data rows

    def test_prefilled_from_samplesheet(self, tmpdir, samplesheet):
        """Generate template with sample_id column pre-filled from samplesheet."""
        out = os.path.join(tmpdir, 'sample_metadata.csv')
        result = subprocess.run(
            ['python', TOOL, 'template', '--samplesheet', samplesheet, '--output-dir', tmpdir],
            capture_output=True, text=True
        )
        assert result.returncode == 0
        with open(out) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            assert len(rows) == 3
            assert rows[0]['sample_id'] == 'ID00001'
            assert rows[1]['sample_id'] == 'ID00002'
            assert rows[0]['organism'] == 'Staphylococcus aureus'

    def test_antibiogram_template(self, tmpdir, samplesheet):
        """Generate antibiogram template alongside metadata template."""
        result = subprocess.run(
            ['python', TOOL, 'template', '--samplesheet', samplesheet,
             '--antibiogram', '--output-dir', tmpdir],
            capture_output=True, text=True
        )
        assert result.returncode == 0
        abg_path = os.path.join(tmpdir, 'antibiogram.csv')
        assert os.path.exists(abg_path)
        with open(abg_path) as f:
            reader = csv.reader(f)
            header = next(reader)
            assert header[0] == 'sample_id'
            assert 'antibiotic' in header
            assert 'resistance_phenotype' in header
