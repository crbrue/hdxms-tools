from __future__ import annotations
from collections import defaultdict
from pathlib import Path
import math
import numpy as np
import pandas as pd
from .io import load_hdx_csv, apply_offsets


def consensus_from_dataframe(df: pd.DataFrame, *, chain: str = "", stat: str = "mean", mincov: int = 1,
                             smooth: int = 0, resi_min: int | None = None, resi_max: int | None = None,
                             value_name: str = "pctD_diff") -> pd.DataFrame:
    if chain:
        df = df[df["chain"] == chain].copy()
    if df.empty:
        raise ValueError("No HDX rows remain after filtering.")
    buckets: dict[int, list[tuple[float, int]]] = defaultdict(list)
    for row in df.itertuples(index=False):
        length = max(1, int(row.end) - int(row.start) + 1)
        for residue in range(int(row.start), int(row.end) + 1):
            buckets[residue].append((float(row.value), length))
    lo = int(resi_min if resi_min is not None else df["start"].min())
    hi = int(resi_max if resi_max is not None else df["end"].max())
    values: dict[int, float] = {}
    coverage: dict[int, int] = {}
    stdev: dict[int, float] = {}
    sem: dict[int, float] = {}
    for residue in range(lo, hi + 1):
        pairs = buckets.get(residue, [])
        nums = np.asarray([v for v, _ in pairs], dtype=float)
        coverage[residue] = len(nums)
        if len(nums) < mincov:
            values[residue] = np.nan; stdev[residue] = np.nan; sem[residue] = np.nan
            continue
        if stat == "median":
            values[residue] = float(np.median(nums))
        elif stat == "len_weighted_mean":
            weights = np.asarray([length for _, length in pairs], dtype=float)
            values[residue] = float(np.average(nums, weights=weights))
        elif stat == "mean":
            values[residue] = float(np.mean(nums))
        else:
            raise ValueError("stat must be mean, median, or len_weighted_mean")
        stdev[residue] = float(np.std(nums, ddof=0)) if len(nums) >= 2 else np.nan
        sem[residue] = stdev[residue] / math.sqrt(len(nums)) if len(nums) >= 2 else np.nan
    if smooth:
        if smooth < 3 or smooth % 2 == 0:
            raise ValueError("smooth must be 0 or an odd integer >= 3")
        series = pd.Series(values, dtype=float)
        smoothed = series.rolling(smooth, center=True, min_periods=1).mean()
        for residue in values:
            values[residue] = float(smoothed.loc[residue]) if coverage[residue] > 0 else np.nan
    return pd.DataFrame({
        "resi": range(lo, hi + 1),
        value_name: [values[r] for r in range(lo, hi + 1)],
        "coverage": [coverage[r] for r in range(lo, hi + 1)],
        "stdev": [stdev[r] for r in range(lo, hi + 1)],
        "sem": [sem[r] for r in range(lo, hi + 1)],
    })


def consensus_from_csv(path: str | Path, *, value_col: str = "pctD_diff", chain: str = "", offset: int = 0,
                       offset_per_chain: str | dict[str, int] | None = None, stat: str = "mean", mincov: int = 1,
                       smooth: int = 0, resi_min: int | None = None, resi_max: int | None = None,
                       output: str | Path | None = None) -> pd.DataFrame:
    df = apply_offsets(load_hdx_csv(path, value_col), offset, offset_per_chain)
    result = consensus_from_dataframe(df, chain=chain, stat=stat, mincov=mincov, smooth=smooth,
                                      resi_min=resi_min, resi_max=resi_max, value_name=value_col)
    if output:
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(output, index=False)
    return result
