from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd


def calculate_back_exchange_stats(combined_csv: str, per_output: str | None = None,
                                  global_output: str | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(combined_csv)
    corrected = [c for c in df.columns if "_emp_corr_" in c and ("Deut %" in c or "% Max" in c)]
    if not corrected:
        raise ValueError("No empirical-corrected deuteration columns were found.")
    back_cols = []
    for corrected_col in corrected:
        raw_col = corrected_col.replace("_emp_corr_", "_")
        if raw_col not in df.columns:
            raise KeyError(f"Expected raw column {raw_col!r} for {corrected_col!r}")
        replicate = corrected_col.split("_emp_corr_", 1)[1]
        out_col = f"back_exch_{replicate}"
        df[out_col] = (df[raw_col] - df[corrected_col]).abs()
        back_cols.append(out_col)
    df["mean_back_exch"] = df[back_cols].mean(axis=1)
    df["std_back_exch"] = df[back_cols].std(axis=1)
    df["median_back_exch"] = df[back_cols].median(axis=1)
    id_cols = [c for c in ["Start", "End", "Sequence", "Charge"] if c in df.columns]
    per_df = df[id_cols + back_cols + ["mean_back_exch", "std_back_exch", "median_back_exch"]]
    rows = []
    for col in back_cols:
        rep = col.removeprefix("back_exch_")
        vals = df[col]
        rows.append({"Condition": rep.split("_noFD")[0] if "_noFD" in rep else rep.split("_")[0],
                     "Replicate": rep, "mean_back_exch": vals.mean(), "std_back_exch": vals.std(),
                     "median_back_exch": vals.median()})
    all_values = df[back_cols].to_numpy(dtype=float).ravel()
    rows.append({"Condition": "All", "Replicate": "All", "mean_back_exch": np.nanmean(all_values),
                 "std_back_exch": np.nanstd(all_values, ddof=1), "median_back_exch": np.nanmedian(all_values)})
    global_df = pd.DataFrame(rows)
    if per_output: Path(per_output).parent.mkdir(parents=True, exist_ok=True); per_df.to_csv(per_output, index=False)
    if global_output: Path(global_output).parent.mkdir(parents=True, exist_ok=True); global_df.to_csv(global_output, index=False)
    return per_df, global_df
