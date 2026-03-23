# Staphit Enhancement Plan: QC & Species Screening

## Motivation
Based on code review (2026-03-23) and PhD workflow needs:

1. **Panaroo risk**: Runs on ALL assemblies, including contaminants
2. **SRA fragility**: Downloads fail silently on network/accession errors
3. **Hardcoded thresholds**: Assembly filter (500KB) not configurable
4. **Limited QC feedback**: Users can't easily diagnose pipeline failures

---

## Proposed Enhancements

### 1. Species Screening with MASH (High Priority)
**Problem**: Panaroo assumes all assemblies are S. aureus. Contaminants cause pangenome failures or garbage output.

**Solution**: Add MASH screening step after assembly, before annotation.

```
ASSEMBLIES → MASH SKETCH → MASH SCREEN (vs S. aureus reference) → FILTER → ANNOTATION
```

**Implementation**:
- Add `S. aureus USA300 (RefSeq NC_007793)` as reference sketch
- Threshold: mash dist < 0.05 (~95% ANI)
- Flag failed samples in summary report
- Allow override: `params.skip_species_check = false`

**Files**: 
- `modules/screen_species.nf` (new)
- `main.nf` (add SCREEN process)

---

### 2. SRA Download Resilience (Medium Priority)
**Problem**: Single SRA download failure kills pipeline. No retry, no fallback.

**Solution**: Add retry wrapper and better error reporting.

**Implementation**:
- Retry up to 3 times with exponential backoff
- Log accession + error to `failed_sra.tsv`
- Continue pipeline with successful samples
- Add `params.max_sra_retries = 3`

**Files**:
- `modules/fetch_sra.nf` (modify)
- `main.nf` (add .catch{} wrapper)

---

### 3. Configurable QC Thresholds (Medium Priority)
**Problem**: 500KB assembly filter is hardcoded. Different species/biases need different thresholds.

**Solution**: Make all QC thresholds parameters.

```groovy
params {
    min_assembly_size = 500000    // bytes
    min_completeness = null        // % (via QUAST)
    max_contamination = null       // % (via CheckM/QUAST)
    max_n50 = null                // or minimum N50
    min_n50 = null
}
```

**Implementation**:
- Add to `params` block in `nextflow.config`
- Update filter logic in `main.nf`
- Document in README.md

**Files**:
- `nextflow.config`
- `main.nf`
- `README.md`

---

### 4. Enhanced Summary Report (Low Priority)
**Problem**: Users must parse multiple CSVs to understand pipeline results.

**Solution**: Generate single HTML summary dashboard.

**Features**:
- Assembly stats (N50, size, completeness)
- Typing results (MLST, Spa, SCCmec, AGR)
- AMR profile (ABRICATE + AMRFinderPlus)
- Species screen results
- Failed samples + reasons
- Downloadable tables

**Implementation**:
- Add `modules/generate_summary.nf` (Python + Jinja2)
- Template: `templates/summary.html.j2`
- Run at pipeline end, output: `results/summary.html`

**Files**:
- `modules/generate_summary.nf` (new)
- `templates/summary.html.j2` (new)
- `main.nf` (call GENERATE_SUMMARY)
- `nextflow.config` (container: `python:3.9-slim`)

---

## Implementation Order

| # | Enhancement | Complexity | Impact | Priority |
|---|------------|------------|--------|----------|
| 1 | Species Screening | Medium | High | P0 |
| 2 | SRA Retry | Low | Medium | P1 |
| 3 | Configurable QC | Low | Medium | P1 |
| 4 | HTML Summary | High | Medium | P2 |

**Sprint 1**: #1 + #2 (MASH screen + SRA retry)
**Sprint 2**: #3 (QC parameters)
**Sprint 3**: #4 (HTML dashboard)

---

## Testing Plan

### Species Screening
```bash
# Test with pure S. aureus (should pass)
nextflow run main.nf -profile docker --input pure_sra_samples.csv

# Test with contaminant mix (should filter)
nextflow run main.nf -profile docker --input mixed_samples.csv --expect_failures 5
```

### SRA Retry
```bash
# Test with invalid accession (should log and continue)
nextflow run main.nf -profile docker --species "Staphylococcus aureus" --search_limit 5

# Check failed_sra.tsv
cat results/failed_sra.tsv
```

### QC Parameters
```bash
# Test strict thresholds
nextflow run main.nf -profile docker --min_assembly_size 1000000 --min_completeness 95

# Verify samples filtered correctly
grep "FILTERED" results/assembly_qc.csv
```

### HTML Summary
```bash
# Run full pipeline
nextflow run main.nf -profile docker

# Open summary
xdg-open results/summary.html
```

---

## Dependencies

**New containers needed**:
- None (using existing staphb/mash + python:3.9-slim)

**Reference genomes**:
- S. aureus USA300 (NC_007793) — already in RefSeq
- Optional: S. epidermidis, S. haemolyticus (for multi-species runs)

**Python packages** (for summary):
- `jinja2` (HTML templating)
- `pandas` (data aggregation)

---

## Success Criteria

- ✅ All assemblies pass species screen before Panaroo
- ✅ Failed SRA downloads logged, pipeline continues
- ✅ QC thresholds configurable via params
- ✅ Single HTML summary generated with all results
- ✅ No breaking changes to existing workflow
- ✅ Documentation updated (README.md, docs/)

---

## Backward Compatibility

All enhancements use new parameters with defaults that maintain existing behavior:
- `params.skip_species_check = false` → enables screen, but user can disable
- `params.max_sra_retries = 3` → added resilience, transparent
- QC params default to null → no change to current filters
- HTML summary → optional addition, doesn't replace existing CSVs

---

**Last Updated**: 2026-03-23
**Branch**: enhancement/qc-and-screening
