# HDX-MS Tools: Rigorous Use, Metadata, and Data Provenance

## 1. Purpose and scientific contract

`hdx-ms-tools` converts HX Examiner CSV exports into reproducible HDX-MS condition summaries, differential tables, consensus residue maps, Illustrator-ready strips, and structure-coloring inputs.

The workflow has two named analysis concepts:

- **`be`**: values normalized against the fully deuterated (FD) control specified in the manifest.
- **`no_be`**: values taken directly from the HX Examiner input exports, with no FD normalization applied.

An empirical D2O-fraction adjustment is **not** a default analysis branch. It is used only when `apply_empirical_correction: true` is explicitly set. The empirical correction must never be inferred from metadata or silently enabled.

The canonical Supplementary Information table must reproduce the original `SI_summary_table_generator.ipynb` output exactly in column names, column order, row order, replicate order, calculations, and rounding behavior. Any additional tables are secondary products.

## 2. Required inputs

Each run requires:

1. One YAML manifest.
2. One fully deuterated control CSV.
3. One or more sample datasets, each containing one or more HX Examiner replicate CSV files.
4. Optional pairwise comparisons.

Replicate order is meaningful. The first file listed is replicate 1, the second is replicate 2, and so on. Matched replicate-difference calculations use this order directly.

### Required HX Examiner columns

Every input must contain:

- `Start`
- `End`
- `Sequence`
- `Charge`
- `# Deut`
- A recognized percent-deuteration column, such as `Deut %`, `% Max D`, `%D`, `Average %D`, `Avg %D`, `Percent D`, `PctD`, or `Pct D`

When present, rows whose `Confidence` value is `Low` are removed before peptide matching.

## 3. Installation

From the package root:

```bash
python -m pip install -e ".[dev,stats]"
```

Confirm the command is available:

```bash
hdxms --help
```

## 4. Recommended project layout

```text
project/
├── exports/
│   ├── A762_FDH2_09222026.csv
│   ├── apo_10s_rep1.csv
│   ├── apo_10s_rep2.csv
│   ├── apo_10s_rep3.csv
│   ├── drugA_10s_rep1.csv
│   ├── drugA_10s_rep2.csv
│   └── drugA_10s_rep3.csv
├── workflow_manifest.yaml
└── hdxms_output/
```

Raw input files should be treated as immutable. Do not overwrite, rename, deduplicate, or edit them after a run has been performed. Any preprocessing should write new files with new names.

## 5. Manifest specification

```yaml
base_dir: .
output_dir: hdxms_output

project:
  project_id: PARP1_HDX_2026_001
  experiment_id: A762_drugA_10s
  operator: "Full name"
  institution: "Institution"
  run_date: "2026-09-22"
  notes: "Optional free-text notes"

sample_d2o_fraction: 0.75
fd_d2o_fraction: 0.75
apply_empirical_correction: false

fully_deuterated: exports/A762_FDH2_09222026.csv

datasets:
  apo_10s:
    title: "Apo, 10 s"
    files:
      - exports/apo_10s_rep1.csv
      - exports/apo_10s_rep2.csv
      - exports/apo_10s_rep3.csv

  drug_a_10s:
    title: "Drug A, 10 s"
    files:
      - exports/drugA_10s_rep1.csv
      - exports/drugA_10s_rep2.csv
      - exports/drugA_10s_rep3.csv

comparisons:
  - name: drug_a_vs_apo_10s
    dataset_a: drug_a_10s
    dataset_b: apo_10s
    title: "PARP1 Δ%HDX — Drug A − Apo, 10 s"
```

### Manifest rules

- Paths may be absolute or relative to `base_dir`.
- Dataset names and comparison names must be unique and stable.
- Replicate files must be listed in the intended biological/technical replicate order.
- `fully_deuterated` must identify the exact FD file used for every `be` normalization in the run.
- A separate manifest and output directory should be used when different datasets require different FD controls.
- `apply_empirical_correction` must remain `false` unless the correction is intentionally requested.
- When empirical correction is enabled, both D2O fractions must be recorded explicitly.

## 6. Running the workflow

Directly:

```bash
hdxms run workflow_manifest.yaml
```

Using the orchestrator:

```bash
bash run_complete_hdx_workflow.sh workflow_manifest.yaml
```

The shell script exits immediately if `hdxms` is not installed or if the workflow fails.

## 7. Provenance outputs

Every run creates a `provenance/` directory. These files are part of the scientific result and should be archived with all figures and tables.

### `provenance/manifest_used.yaml`

An exact copy of the manifest used for the run. This prevents later edits to the working manifest from obscuring the original configuration.

### `provenance/raw_input_inventory.csv`

One row per raw input file, containing:

- input role (`sample` or `fully_deuterated`)
- dataset name
- replicate number
- exact input filename
- resolved absolute input path
- SHA-256 checksum
- file size in bytes
- filesystem modification timestamp in UTC
- exact FD-control filename used for normalization
- resolved FD-control path
- FD-control SHA-256 checksum
- sample D2O fraction
- FD-control D2O fraction

The SHA-256 checksum is the authoritative identity of a file. Two files with the same name but different checksums are different inputs.

### `provenance/peptide_attrition_by_input.csv`

One row per input file. It reports both unique peptide sequences and peptide/charge pairs:

- count before common-peptide filtering
- count after common-peptide filtering
- count lost during common-peptide filtering
- count retained after `be` normalization
- count lost during `be` normalization

A peptide/charge pair is defined by `(Sequence, Charge)`. This is the unit used for matching because the same sequence at different charge states is analytically distinct.

### `provenance/peptide_attrition_by_dataset.csv`

One row per dataset, containing:

- replicate count
- union of peptide/charge pairs observed across replicates
- common peptide/charge pairs shared within that dataset before global filtering
- pairs retained in the global common set
- pairs lost when enforcing the global common set
- common pairs still valid after FD normalization
- pairs lost during FD normalization
- FD-control filename and checksum

### `provenance/software_environment.yaml`

Records the Python, operating-system, pandas, NumPy, and PyYAML versions used for the run.

### `run_metadata.yaml`

The run-level index. It records:

- UTC run timestamp
- manifest path and checksum
- output directory
- primary normalization mode
- exact FD-control name, path, and checksum
- D2O fractions
- empirical-correction status
- global common peptide/charge count
- global common unique-sequence count
- every dataset and its ordered replicate filenames
- the provenance files generated

## 8. Peptide attrition definitions

Peptide loss must be reported at two levels.

### A. Loss caused by common-peptide matching

For an input file:

```text
lost_common = pairs_before_filter − pairs_in_global_common_set
```

The global common set is the intersection of `(Sequence, Charge)` pairs across every sample replicate and the FD control included in the run.

This is intentionally strict. A pair absent from any included file is excluded from all downstream comparisons.

### B. Loss caused by FD normalization

For a filtered input file:

```text
lost_be = common_pairs − pairs_with_finite_normalized_values
```

A pair is counted as lost during FD normalization if the resulting normalized `# Deut` or percent-deuteration value is missing or non-finite, for example because the FD denominator is missing, zero, or invalid.

The attrition report distinguishes these losses. A peptide removed because it is not common is not counted again as an FD-normalization loss.

## 9. Required review before accepting a run

A run should not be accepted until the following checks are completed:

1. Verify every expected raw filename appears in `raw_input_inventory.csv`.
2. Verify each sample row points to the intended FD-control filename and checksum.
3. Verify replicate numbers match the order in the manifest.
4. Review peptide losses in `peptide_attrition_by_input.csv`.
5. Review dataset-level losses in `peptide_attrition_by_dataset.csv`.
6. Investigate any nonzero `pairs_lost_during_be_normalization` value.
7. Confirm the SI table has the exact legacy column order and replicate order.
8. Confirm comparison direction: `pctD_diff = dataset_a − dataset_b`.
9. Archive the manifest, provenance directory, SI tables, comparison tables, and figures together.

## 10. Reproducibility and reruns

Never overwrite an accepted output directory. Use a unique run directory, for example:

```yaml
output_dir: outputs/PARP1_HDX_2026_001_run01
```

For a revised run:

```yaml
output_dir: outputs/PARP1_HDX_2026_001_run02
```

Differences between runs should be evaluated using:

- manifest checksums
- raw-input checksums
- FD-control checksum
- software environment
- peptide attrition reports
- final table checksums, if externally archived

## 11. Interpretation of `be`, `no_be`, and empirical correction

### `be`

Uses the FD control specified in `fully_deuterated`. Every `be` table and figure must be traceable to that exact FD filename and checksum.

### `no_be`

Uses the values supplied in each HX Examiner input without applying FD normalization. The FD file may still participate in global common-peptide matching only when the workflow is configured to use one common set across both branches; this behavior must be recorded in metadata.

### Empirical correction

Applied only when explicitly enabled. It is not synonymous with `be`, and it must not replace or silently modify the standard `be` branch. The manifest must record both D2O fractions and the enabled flag.

## 12. Minimum archive for publication or handoff

Archive the following together:

```text
workflow_manifest.yaml
provenance/manifest_used.yaml
provenance/raw_input_inventory.csv
provenance/peptide_attrition_by_input.csv
provenance/peptide_attrition_by_dataset.csv
provenance/software_environment.yaml
run_metadata.yaml
all canonical SI tables
all differential tables
all consensus tables
all final figures and SVG files
package version or source-code commit
```

Without this set, the analysis should be considered incomplete from a provenance standpoint.

## 11. Canonical dual SI-table requirement

Every dataset must produce two structurally identical canonical SI tables:

```text
datasets/<dataset>/be/SI_summary_table.csv
datasets/<dataset>/no_be/SI_summary_table.csv
```

The `be` table applies the selected FD/back-exchange normalization. The `no_be` table contains the raw HX Examiner uptake values after only documented filtering and peptide matching. No QC or provenance columns may be added to either canonical SI table.

## 12. Back-exchange audit data

Per-peptide back-exchange values and the median applied to the BE branch must be exported separately under:

```text
datasets/<dataset>/back_exchange/
```

See `OUTPUT_CONTRACT.md` for the required columns and `TESTING.md` for manual validation steps.

## 13. Peptide-level QC

Peptide QC is generated after the scientific output tables and is diagnostic only:

```text
analysis_report/peptide_quality_scores.csv
analysis_report/red_flagged_peptides.csv
```

The red-flag table is the subset of peptides categorized as `Poor`. QC scoring must never silently remove peptides from the SI tables or comparisons. See `PEPTIDE_QC_SPEC.md`.

## 14. Analysis report

Every run must produce:

```text
analysis_report/analysis_report.yaml
```

This file summarizes data quality, peptide attrition, FD behavior, back exchange, peptide-level QC, warnings, and the generated-output inventory. See `ANALYSIS_REPORT_SPEC.md`.
