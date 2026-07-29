# Peptide-Level Quality Control Specification

## 1. Purpose

Peptide QC provides a transparent diagnostic assessment for every peptide/charge pair. It does not automatically exclude data.

## 2. Required outputs

```text
analysis_report/peptide_quality_scores.csv
analysis_report/red_flagged_peptides.csv
```

## 3. Required columns

The complete table should include:

- `Dataset`
- `Start`
- `End`
- `Sequence`
- `Charge`
- `FD_source_file`
- `raw_FD_uptake`
- `theoretical_maxD`
- `FD_fraction_of_theoretical_max`
- `per_peptide_back_exchange_fraction`
- `median_back_exchange_fraction_applied`
- `maximum_no_be_uptake`
- `maximum_be_uptake`
- `maximum_corrected_occupancy_fraction`
- `mean_replicate_SD_deuterons`
- `maximum_replicate_SD_deuterons`
- `missing_replicate_count`
- `missing_timepoint_count`
- `QC_score`
- `QC_category`
- `QC_flags`
- `Recommendation`

## 4. Core flags

### Corrected occupancy

- `FD_OVERCORRECTION_GT_100`: maximum corrected occupancy is greater than 1.00
- `FD_OVERCORRECTION_GT_110`: maximum corrected occupancy is greater than 1.10

The severe flag supersedes the warning flag for scoring. Do not double-penalize the same event.

### Replicate variability

- `HIGH_REPLICATE_SD`
- `VERY_HIGH_REPLICATE_SD`

Thresholds are configurable in the manifest and are expressed in absolute deuterons unless a separate percent-based metric is explicitly added.

### FD validity

- `MISSING_FD_MATCH`
- `DUPLICATE_FD_MATCH`
- `NONFINITE_FD_VALUE`
- `FD_UPTAKE_GT_THEORETICAL_MAX`
- `FD_UPTAKE_LE_ZERO`

### Sample completeness

- `MISSING_REPLICATE`
- `MISSING_TIMEPOINT`
- `NONFINITE_UPTAKE`

## 5. Default scoring model

Start every peptide at 100 points. Apply only the most severe penalty within each nested flag family.

Suggested defaults:

| Condition | Penalty |
|---|---:|
| Corrected occupancy >100% | 10 |
| Corrected occupancy >110% | 25 |
| High replicate SD | 10 |
| Very high replicate SD | 25 |
| One missing replicate | 15 |
| Multiple missing replicates | 30 |
| Missing FD match | 30 |
| Invalid FD value | 20 |

Clamp scores to the range 0–100.

Suggested categories:

- `Excellent`: 90–100
- `Acceptable`: 75–89
- `Review`: 50–74
- `Poor`: 0–49

## 6. Red-flag export

`red_flagged_peptides.csv` contains only rows where:

```text
QC_category == Poor
```

It should retain all columns from the master QC table so that it remains a fully traceable filtered view rather than a separately maintained dataset.

## 7. Recommendations

Recommendations should be derived from flags. Examples:

- `FD_OVERCORRECTION_GT_110`: review FD assignment, theoretical maximum, and FD intensity quality
- `VERY_HIGH_REPLICATE_SD`: inspect replicate chromatograms and uptake extraction
- `MISSING_FD_MATCH`: verify FD export and peptide identity matching
- `NONFINITE_UPTAKE`: inspect raw HX Examiner export for missing or malformed values
- `MISSING_REPLICATE`: confirm replicate export completeness

## 8. Manifest configuration

```yaml
quality_control:
  enabled: true

  corrected_occupancy:
    warning: 1.00
    severe: 1.10

  replicate_sd_deuterons:
    warning: 0.50
    severe: 1.00

  scoring:
    corrected_occupancy_warning_penalty: 10
    corrected_occupancy_severe_penalty: 25
    high_sd_penalty: 10
    severe_sd_penalty: 25
    one_missing_replicate_penalty: 15
    multiple_missing_replicates_penalty: 30
    missing_fd_penalty: 30
    invalid_fd_penalty: 20

  categories:
    excellent_min: 90
    acceptable_min: 75
    review_min: 50
```
