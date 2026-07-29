from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

KEYS = ["Start", "End", "Sequence", "Charge"]


def theoretical_max_d(sequence: str) -> int:
    """Notebook-compatible exchangeable-amide estimate."""
    seq = str(sequence)
    return max(len(seq) - 1 - seq.count("P"), 0)


def _pct_col(df: pd.DataFrame) -> str:
    candidates = ["Deut %", "% Max D", "%D", "Average %D", "Avg %D", "Percent D", "PctD", "Pct D"]
    for col in candidates:
        if col in df.columns:
            return col
    for col in df.columns:
        text = str(col).lower()
        if "%d" in text or "percent" in text or ("deut" in text and "%" in str(col)):
            return col
    raise KeyError("Could not identify a percent-deuteration column")


def build_back_exchange_table(fd: pd.DataFrame, fd_fraction: float) -> pd.DataFrame:
    """Create an auditable, per-peptide FD/back-exchange table.

    Back exchange is reported as 1 - measured_FD_D / expected_FD_D, where
    expected_FD_D = theoretical_max_D * fd_fraction.
    """
    out = fd[KEYS].copy()
    out["theoretical_maxD"] = out["Sequence"].map(theoretical_max_d)
    out["fd_d2o_fraction"] = float(fd_fraction)
    out["expected_fd_D"] = out["theoretical_maxD"] * float(fd_fraction)
    out["raw_fd_uptake_D"] = pd.to_numeric(fd["# Deut"], errors="coerce")
    out["raw_fd_percent_column"] = pd.to_numeric(fd[_pct_col(fd)], errors="coerce")
    denom = out["expected_fd_D"].replace(0, np.nan)
    out["fd_occupancy_fraction"] = out["raw_fd_uptake_D"] / denom
    out["peptide_back_exchange_fraction"] = 1.0 - out["fd_occupancy_fraction"]
    valid = out["peptide_back_exchange_fraction"].replace([np.inf, -np.inf], np.nan)
    median = float(valid.median()) if valid.notna().any() else float("nan")
    out["median_back_exchange_fraction"] = median
    out["included_in_median"] = valid.notna()
    return out.sort_values(["Start", "End", "Charge"]).reset_index(drop=True)


def build_peptide_quality_table(
    dataset: str,
    raw_frames: list[pd.DataFrame],
    be_frames: list[pd.DataFrame],
    fd_table: pd.DataFrame,
    qc_cfg: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """Score peptide/charge pairs diagnostically; never exclude them."""
    qc_cfg = qc_cfg or {}
    occ_cfg = qc_cfg.get("corrected_occupancy", {})
    sd_cfg = qc_cfg.get("replicate_sd_deuterons", {})
    scoring = qc_cfg.get("scoring", {})
    categories = qc_cfg.get("categories", {})

    occ_warn = float(occ_cfg.get("warning", 1.00))
    occ_severe = float(occ_cfg.get("severe", 1.10))
    sd_warn = float(sd_cfg.get("warning", 0.50))
    sd_severe = float(sd_cfg.get("severe", 1.00))

    fd_lookup = fd_table.set_index(["Sequence", "Charge"])
    base = raw_frames[0][KEYS].copy().reset_index(drop=True)
    base["Dataset"] = dataset
    base["theoretical_maxD"] = base["Sequence"].map(theoretical_max_d)

    raw_values = np.column_stack([pd.to_numeric(f["# Deut"], errors="coerce").to_numpy() for f in raw_frames])
    be_values = np.column_stack([pd.to_numeric(f["# Deut"], errors="coerce").to_numpy() for f in be_frames])
    base["replicate_count_expected"] = len(raw_frames)
    base["replicate_count_observed"] = np.isfinite(raw_values).sum(axis=1)
    base["mean_raw_uptake_D"] = np.nanmean(raw_values, axis=1)
    base["max_raw_uptake_D"] = np.nanmax(raw_values, axis=1)
    base["replicate_sd_D"] = np.nanstd(raw_values, axis=1, ddof=1) if len(raw_frames) > 1 else 0.0
    base["max_be_corrected_uptake"] = np.nanmax(be_values, axis=1)
    denom = base["theoretical_maxD"].replace(0, np.nan)
    base["max_corrected_occupancy"] = base["max_be_corrected_uptake"] / denom

    fd_uptake, fd_occ, be_frac = [], [], []
    for seq, charge in zip(base["Sequence"], base["Charge"]):
        key = (str(seq), charge)
        if key in fd_lookup.index:
            row = fd_lookup.loc[key]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            fd_uptake.append(row.get("raw_fd_uptake_D", np.nan))
            fd_occ.append(row.get("fd_occupancy_fraction", np.nan))
            be_frac.append(row.get("peptide_back_exchange_fraction", np.nan))
        else:
            fd_uptake.append(np.nan); fd_occ.append(np.nan); be_frac.append(np.nan)
    base["raw_fd_uptake_D"] = fd_uptake
    base["fd_occupancy_fraction"] = fd_occ
    base["peptide_back_exchange_fraction"] = be_frac

    scores, labels, flags_col, recommendations = [], [], [], []
    for row in base.itertuples(index=False):
        score = 100
        flags: list[str] = []
        recs: list[str] = []
        occ = getattr(row, "max_corrected_occupancy")
        sd = getattr(row, "replicate_sd_D")
        missing = int(getattr(row, "replicate_count_expected") - getattr(row, "replicate_count_observed"))
        fd_value = getattr(row, "raw_fd_uptake_D")

        if not np.isfinite(fd_value):
            flags.append("MISSING_FD_MATCH")
            recs.append("Verify FD export and peptide/charge matching")
            score -= int(scoring.get("missing_fd_penalty", 30))
        elif fd_value <= 0:
            flags.append("INVALID_FD_VALUE")
            recs.append("Inspect the FD uptake value")
            score -= int(scoring.get("invalid_fd_penalty", 20))

        if np.isfinite(occ) and occ > occ_severe:
            flags.append("FD_OVERCORRECTION_GT_110")
            recs.append("Review FD assignment and theoretical maximum deuterons")
            score -= int(scoring.get("corrected_occupancy_severe_penalty", 30))
        elif np.isfinite(occ) and occ > occ_warn:
            flags.append("FD_OVERCORRECTION_GT_100")
            recs.append("Review corrected occupancy")
            score -= int(scoring.get("corrected_occupancy_warning_penalty", 10))

        if np.isfinite(sd) and sd > sd_severe:
            flags.append("VERY_HIGH_REPLICATE_SD")
            recs.append("Inspect raw replicate uptake values for outliers")
            score -= int(scoring.get("severe_sd_penalty", 30))
        elif np.isfinite(sd) and sd > sd_warn:
            flags.append("HIGH_REPLICATE_SD")
            recs.append("Inspect replicate variability")
            score -= int(scoring.get("high_sd_penalty", 10))

        if missing > 1:
            flags.append("MULTIPLE_MISSING_REPLICATES")
            recs.append("Confirm replicate export completeness")
            score -= int(scoring.get("multiple_missing_replicates_penalty", 30))
        elif missing == 1:
            flags.append("MISSING_REPLICATE")
            recs.append("Confirm replicate export completeness")
            score -= int(scoring.get("one_missing_replicate_penalty", 15))

        score = max(0, score)
        if score >= int(categories.get("excellent_min", 90)):
            label = "Excellent"
        elif score >= int(categories.get("acceptable_min", 75)):
            label = "Acceptable"
        elif score >= int(categories.get("review_min", 50)):
            label = "Review"
        else:
            label = "Poor"
        scores.append(score)
        labels.append(label)
        flags_col.append(";".join(flags) if flags else "NONE")
        recommendations.append("; ".join(dict.fromkeys(recs)) if recs else "No action required")

    base["QC Score"] = scores
    base["QC Category"] = labels
    base["QC Flags"] = flags_col
    base["Recommendation"] = recommendations
    return base


def write_quality_outputs(table: pd.DataFrame, output_dir: str | Path) -> tuple[Path, Path]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    all_path = output_dir / "peptide_quality_scores.csv"
    red_path = output_dir / "flagged_peptides.csv"
    table.to_csv(all_path, index=False)
    table.loc[table["QC Category"] == "Poor"].to_csv(red_path, index=False)
    return all_path, red_path


def quality_summary(table: pd.DataFrame) -> dict[str, Any]:
    counts = table["QC Category"].value_counts().to_dict()
    return {
        "total_peptide_charge_pairs": int(len(table)),
        "excellent": int(counts.get("Excellent", 0)),
        "acceptable": int(counts.get("Acceptable", 0)),
        "review": int(counts.get("Review", 0)),
        "poor": int(counts.get("Poor", 0)),
        "overcorrected_gt_100": int(table["QC Flags"].str.contains("GT_100", regex=False).sum()),
        "overcorrected_gt_110": int(table["QC Flags"].str.contains("GT_110", regex=False).sum()),
        "high_sd": int(table["QC Flags"].str.contains("HIGH_REPLICATE_SD", regex=False).sum()),
        "very_high_sd": int(table["QC Flags"].str.contains("VERY_HIGH_REPLICATE_SD", regex=False).sum()),
    }
