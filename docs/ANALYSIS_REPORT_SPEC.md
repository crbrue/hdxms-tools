# Analysis Report Specification

## 1. Required file

```text
analysis_report/analysis_report.yaml
```

The report is the primary run-level QC and provenance summary.

## 2. Required sections

### `run`

- UTC timestamp
- project ID
- experiment ID
- package version
- git commit, when available
- manifest path and SHA-256
- output directory
- Python version
- operating system

### `inputs`

For each raw file:

- role
- dataset
- replicate number
- exact filename
- resolved path
- SHA-256
- file size
- UTC modification timestamp

### `fd_control`

- exact filename
- resolved path
- SHA-256
- sample D2O fraction
- FD D2O fraction
- empirical correction enabled or disabled

### `datasets`

For each dataset:

- title
- replicate count
- raw peptide/charge-pair union
- within-dataset common peptide/charge pairs
- globally retained peptide/charge pairs
- pairs lost during common filtering
- pairs retained after FD normalization
- pairs lost during FD normalization
- unique peptide sequences
- mean and median replicate SD

### `back_exchange`

For each dataset:

- number of valid FD-matched peptide/charge pairs
- number excluded from the median calculation
- mean back exchange
- median back exchange
- standard deviation
- IQR
- minimum
- maximum
- exact median applied

### `peptide_quality`

For each dataset and globally:

- total evaluated
- Excellent count
- Acceptable count
- Review count
- Poor count
- corrected occupancy >100% count
- corrected occupancy >110% count
- high replicate SD count
- very high replicate SD count
- missing FD match count

### `comparisons`

For each comparison:

- comparison direction
- final peptide/charge-pair count
- significant peptide count, when significance criteria are configured
- largest positive difference
- largest negative difference

### `warnings`

Warnings must be explicit and machine-readable. Each warning should include:

- code
- severity
- dataset or comparison
- message
- relevant count or value
- recommended review action

### `outputs`

List every expected output with:

- path
- status: `created`, `skipped`, or `failed`
- reason when skipped or failed

## 3. Quality interpretation

An overall run category may be reported, but only if the criteria are explicit. The report must always expose the underlying metrics and warnings. A single score must never replace the detailed QC evidence.
