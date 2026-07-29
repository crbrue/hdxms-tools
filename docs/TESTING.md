# Testing Guide

## 1. Install

```bash
cd hdx-ms-tools
python -m pip install -e ".[dev,stats]"
```

## 2. Run unit tests

```bash
pytest -q
```

## 3. Run a manifest

```bash
hdxms run examples/workflow_manifest.yaml
```

## 4. Acceptance checklist

### Canonical SI tables

Confirm every dataset contains:

```text
datasets/<dataset>/be/SI_summary_table.csv
datasets/<dataset>/no_be/SI_summary_table.csv
```

Confirm the BE and NO_BE files have identical:

- headers
- column order
- row count
- peptide order
- replicate order

Confirm NO_BE values match the raw HX Examiner values after only the documented row filtering and common-set matching.

Confirm BE values reproduce the original notebook output.

### Back exchange

Confirm every dataset contains:

```text
datasets/<dataset>/back_exchange/peptide_back_exchange.csv
datasets/<dataset>/back_exchange/back_exchange_summary.yaml
```

Manually verify several peptides:

1. theoretical maximum deuterons
2. raw FD uptake
3. FD fraction of theoretical maximum
4. per-peptide back-exchange fraction
5. median back-exchange fraction applied

Confirm the median applied in the BE workflow exactly matches the value reported in the summary.

### Peptide QC

Confirm:

```text
analysis_report/peptide_quality_scores.csv
analysis_report/red_flagged_peptides.csv
```

Create or identify test peptides that trigger:

- corrected occupancy >100%
- corrected occupancy >110%
- high replicate SD
- very high replicate SD
- missing FD match

Confirm `red_flagged_peptides.csv` is exactly the subset where `QC_category == Poor`.

Confirm no peptide is automatically removed from the SI table solely because of its QC score.

### Provenance

Confirm every input file appears in `provenance/raw_input_inventory.csv` with the correct SHA-256.

Confirm the exact FD filename and checksum are recorded for each dataset.

Confirm peptide losses are separated into:

- loss during common filtering
- loss during FD normalization

### Analysis report

Confirm `analysis_report/analysis_report.yaml` contains:

- run identity
- software environment
- raw inputs
- FD control
- peptide attrition
- back-exchange statistics
- QC category counts
- warnings
- output inventory

## 5. Current implementation note

At the time this document was written, the package's legacy workflow code may not yet materialize every output in this specification. Treat missing BE/NO_BE parallel directories, back-exchange audit tables, peptide QC tables, or the full analysis report as implementation failures to be addressed, not as optional outputs.
