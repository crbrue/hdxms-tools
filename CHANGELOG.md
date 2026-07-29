# Changelog

## 1.2.0 — Pyteomics Mass Restoration and Modified Peptides

- Restored Pyteomics as a required core dependency for SI-table monoisotopic-mass calculations.
- Preserved the legacy SI-table convention: neutral monoisotopic mass plus `charge × proton mass`.
- Added explicit validation for higher charge states to prevent accidental m/z reporting.
- Added modification-aware mass calculation when `Peptide ID` contains valid ProForma notation.
- Added strict Sequence/Peptide-ID consistency checks; modified identifiers are never silently ignored.
- Added mass regression tests for unmodified peptides, charge states 1–5, and an oxidized peptide.
- Added `environment.yml` and `requirements.txt` with Pyteomics as a core dependency.

## 1.1.0 — Provenance and Quality Control

- Added parallel `datasets/<dataset>/no_be/SI_summary_table.csv` and `be/SI_summary_table.csv` outputs with identical schemas.
- Added per-peptide FD/back-exchange audit export.
- Added peptide-level diagnostic QC scoring and a red-only `flagged_peptides.csv` subset.
- Added machine-readable `analysis_report/analysis_report.yaml`.
- Expanded raw-input fingerprints, peptide attrition, software-environment, and output-inventory provenance.
- Preserved legacy workflow products and comparison plots.

## 1.1.1

### Added
- `scripts/list_raw_csvs.py` for generating a copy/paste-ready Python list of raw CSV filenames.
- Optional conservative manifest-skeleton inference with explicit reporting of ambiguous files and multiple FD candidates.
