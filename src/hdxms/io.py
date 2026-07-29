from __future__ import annotations
from pathlib import Path
import pandas as pd


def load_hdx_csv(path: str | Path, value_col: str = "pctD_diff") -> pd.DataFrame:
    """Normalize peptide- or residue-level HDX CSV to chain/start/end/value."""
    df = pd.read_csv(path)
    cols = {str(c).lower(): c for c in df.columns}
    vkey = cols.get(value_col.lower())
    if vkey is None:
        raise ValueError(f"Value column {value_col!r} not found. Columns: {list(df.columns)}")
    if {"resi_start", "resi_end"}.issubset(cols):
        skey, ekey = cols["resi_start"], cols["resi_end"]
    elif {"start", "end"}.issubset(cols):
        skey, ekey = cols["start"], cols["end"]
    elif "resi" in cols:
        skey = ekey = cols["resi"]
    else:
        raise ValueError("CSV requires Start/End, resi_start/resi_end, or resi columns.")
    chain = df[cols["chain"]].fillna("").astype(str).str.strip() if "chain" in cols else pd.Series("", index=df.index)
    out = pd.DataFrame({
        "chain": chain,
        "start": pd.to_numeric(df[skey], errors="raise").astype(int),
        "end": pd.to_numeric(df[ekey], errors="raise").astype(int),
        "value": pd.to_numeric(df[vkey], errors="raise").astype(float),
    })
    swap = out["end"] < out["start"]
    out.loc[swap, ["start", "end"]] = out.loc[swap, ["end", "start"]].to_numpy()
    return out


def parse_chain_offsets(spec: str | dict[str, int] | None) -> dict[str, int]:
    if spec is None:
        return {}
    if isinstance(spec, dict):
        return {str(k): int(v) for k, v in spec.items()}
    result: dict[str, int] = {}
    for item in str(spec).split(","):
        if item.strip():
            chain, value = item.split(":", 1)
            result[chain.strip()] = int(value)
    return result


def apply_offsets(df: pd.DataFrame, global_offset: int = 0, per_chain: str | dict[str, int] | None = None) -> pd.DataFrame:
    result = df.copy()
    offsets = parse_chain_offsets(per_chain)
    row_offsets = result["chain"].map(offsets).fillna(global_offset).astype(int)
    result["start"] += row_offsets
    result["end"] += row_offsets
    result = result[result["end"] >= 1].copy()
    result["start"] = result["start"].clip(lower=1)
    return result
