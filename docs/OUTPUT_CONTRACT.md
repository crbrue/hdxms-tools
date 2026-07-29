# Output Contract

This document defines the required output structure for `hdx-ms-tools`.

## 1. Canonical SI tables

The original SI summary table is the primary scientific output and must remain unchanged in structure.

For every dataset, generate two parallel copies:

```text
datasets/<dataset>/
├── be/
│   └── SI_summary_table.csv
└── no_be/
    └── SI_summary_table.csv
```

The two SI tables must have identical:

- column names
- column order
- row order
- peptide ordering
- replicate ordering
- calculations other than normalization
- rounding behavior

Only the numerical uptake values differ:

- `be`: standard fully deuterated-control normalization
- `no_be`: raw HX Examiner uptake values with no FD normalization

No QC, provenance, back-exchange, or warning columns may be added to these canonical SI tables.

## 2. Back-exchange exports

For every dataset, write:

```text
datasets/<dataset>/back_exchange/
├── peptide_back_exchange.csv
├── back_exchange_distribution.csv
└── back_exchange_summary.yaml
```

`peptide_back_exchange.csv` must contain one row per peptide/charge pair and include at least:

- `Start`
- `End`
- `Sequence`
- `Charge`
- `FD_source_file`
- `FD_source_sha256`
- `raw_FD_uptake`
- `theoretical_maxD`
- `FD_fraction_of_theoretical_max`
- `per_peptide_back_exchange_fraction`
- `median_back_exchange_fraction_applied`
- `included_in_median`
- `exclusion_reason`

The per-peptide back-exchange value is diagnostic. The canonical BE branch applies the global median back-exchange correction selected by the workflow, not a different correction to each peptide unless explicitly configured in a future method.

`back_exchange_summary.yaml` must record:

- total peptide/charge pairs
- valid peptide/charge pairs
- excluded peptide/charge pairs
- mean back exchange
- median back exchange
- standard deviation
- interquartile range
- minimum
- maximum
- exact median value applied
- FD filename and SHA-256

## 3. Peptide-level QC

Write:

```text
analysis_report/
├── peptide_quality_scores.csv
└── red_flagged_peptides.csv
```

`peptide_quality_scores.csv` is the complete peptide/charge-level QC table.

`red_flagged_peptides.csv` is an automatically generated subset containing only peptides classified as `Poor`.

The QC output is diagnostic only. Peptides must not be silently removed from SI tables or comparisons because of their QC score.

## 4. Analysis report

Write machine-readable run-level QC and provenance to:

```text
analysis_report/analysis_report.yaml
```

The report must summarize:

- run identity and software environment
- exact raw inputs and checksums
- exact FD control and checksum
- dataset replicate counts
- peptide attrition
- back-exchange distribution
- peptide QC category counts
- numbers of peptides above 100% and 110% corrected occupancy
- high-SD peptide counts
- comparison-level peptide counts
- output inventory
- warnings

## 5. Provenance

Every run must include:

```text
provenance/
├── manifest_used.yaml
├── raw_input_inventory.csv
├── peptide_attrition_by_input.csv
├── peptide_attrition_by_dataset.csv
└── software_environment.yaml
```

All raw inputs must be identified by exact filename, resolved path, file size, modification timestamp, and SHA-256 checksum.

## 6. Empirical D2O correction

Empirical D2O correction is never inferred and is never a default branch. It is performed only when explicitly enabled in the manifest.

It must not replace or redefine the canonical `be` and `no_be` outputs without an explicit manifest request and clear metadata indicating the applied method.
