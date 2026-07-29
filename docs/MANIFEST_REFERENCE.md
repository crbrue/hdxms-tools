# Manifest Reference

## Minimal example

```yaml
base_dir: .
output_dir: hdxms_output

project:
  project_id: PARP1_HDX_2026_001
  experiment_id: A762_drugA
  operator: "Full name"
  institution: "Institution"
  run_date: "2026-09-22"
  notes: "Optional notes"

fully_deuterated: exports/A762_FD.csv
sample_d2o_fraction: 0.75
fd_d2o_fraction: 0.75
apply_empirical_correction: false

datasets:
  apo:
    title: "Apo"
    files:
      - exports/apo_rep1.csv
      - exports/apo_rep2.csv
      - exports/apo_rep3.csv

  drug_a:
    title: "Drug A"
    files:
      - exports/drug_a_rep1.csv
      - exports/drug_a_rep2.csv
      - exports/drug_a_rep3.csv

comparisons:
  - name: drug_a_vs_apo
    dataset_a: drug_a
    dataset_b: apo
    title: "Drug A − Apo"

quality_control:
  enabled: true
  corrected_occupancy:
    warning: 1.00
    severe: 1.10
  replicate_sd_deuterons:
    warning: 0.50
    severe: 1.00
```

## Required rules

- Relative paths are resolved from `base_dir`.
- Replicate order is the listed order.
- `dataset_a - dataset_b` defines comparison direction.
- The FD file must be explicitly named.
- Empirical D2O correction remains disabled unless explicitly requested.
- QC thresholds may be customized, but all applied values must be written to `analysis_report.yaml`.
- The workflow always produces both canonical `be` and `no_be` SI tables.
