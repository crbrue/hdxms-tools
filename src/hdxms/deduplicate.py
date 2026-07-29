from __future__ import annotations
from pathlib import Path
import shutil
import pandas as pd


def remove_duplicates(path: str | Path, *, key_cols=("Start", "End", "Sequence", "Charge"),
                      search_rt_col: str = "Search RT", backup: bool = True,
                      backup_suffix: str = ".bak", inplace: bool = True,
                      output: str | Path | None = None) -> pd.DataFrame:
    path = Path(path)
    df = pd.read_csv(path)
    missing = [c for c in (*key_cols, search_rt_col) if c not in df.columns]
    if missing:
        raise KeyError(f"{path}: missing columns {missing}")
    if "Sequence" in key_cols:
        df["Sequence"] = df["Sequence"].astype(str).str.strip()
    for col in set(key_cols).intersection({"Start", "End", "Charge"}):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df[search_rt_col] = pd.to_numeric(df[search_rt_col], errors="coerce")
    ordered = df.assign(_rt_missing=df[search_rt_col].isna()).sort_values(
        list(key_cols) + ["_rt_missing", search_rt_col], kind="stable", na_position="last"
    )
    result = ordered.drop_duplicates(list(key_cols), keep="first").drop(columns="_rt_missing")
    result = result.sort_values(list(key_cols) + [search_rt_col], kind="stable", na_position="last")
    target = path if inplace else Path(output or path.with_name(path.stem + "_deduplicated.csv"))
    if inplace and backup:
        backup_path = Path(str(path) + backup_suffix)
        if not backup_path.exists():
            shutil.copy2(path, backup_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(target) + ".tmp")
    result.to_csv(tmp, index=False)
    tmp.replace(target)
    return result
