from __future__ import annotations

"""Illustrator-safe consensus-strip exporters.

These functions intentionally preserve the manual SVG implementation from
``consensus_plots_forIllustrator_02042026.ipynb``: exact physical dimensions,
one rectangle per residue, a white background for uncovered regions, and
``preserveAspectRatio=\"none\"`` so the strip imports into Illustrator at the
requested dimensions without Matplotlib padding or clipping.
"""

from pathlib import Path
import pandas as pd

from .colors import hex_color


def _load_consensus(csv_path: str, resi_col: str, val_col: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path).copy()
    if resi_col not in df.columns:
        raise ValueError(f"Residue column '{resi_col}' not found. Columns: {list(df.columns)}")
    if val_col not in df.columns:
        raise ValueError(f"Value column '{val_col}' not found. Columns: {list(df.columns)}")
    df[resi_col] = pd.to_numeric(df[resi_col], errors="raise").astype(int)
    df[val_col] = pd.to_numeric(df[val_col], errors="coerce")
    return df.sort_values(resi_col)


def _write_manual_svg(
    *,
    residues,
    values,
    origin: int,
    span: int,
    out_svg: str,
    width_in: float,
    height_in: float,
) -> None:
    if span < 1:
        raise ValueError("SVG residue span must be at least one residue")

    svg = [
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{width_in}in" height="{height_in}in" '
            f'viewBox="0 0 {span} 1" preserveAspectRatio="none" version="1.1">'
        ),
        f'<rect x="0" y="0" width="{span}" height="1" fill="#ffffff" stroke="none"/>',
    ]

    for resi, value in zip(residues, values):
        x = int(resi) - int(origin)
        if 0 <= x < span:
            svg.append(
                f'<rect x="{x}" y="0" width="1" height="1" '
                f'fill="{hex_color(value, missing="#ffffff")}" stroke="none"/>'
            )

    svg.append("</svg>")
    path = Path(out_svg)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(svg), encoding="utf-8")


def export_consensus_full_manual_fixed_size(
    csv_path: str,
    out_svg: str,
    width_in: float = 7.4558,
    height_in: float = 0.0751,
    resi_col: str = "resi",
    val_col: str = "pctD_diff",
) -> None:
    """Export the full Illustrator-ready consensus strip from the notebook."""
    df = _load_consensus(csv_path, resi_col, val_col)
    if df.empty:
        raise ValueError("Consensus CSV contains no rows")

    xs = df[resi_col].to_numpy(dtype=int)
    vals = df[val_col].to_numpy()
    xmin = int(xs.min())
    xmax = int(xs.max())
    span = xmax - xmin + 1

    _write_manual_svg(
        residues=xs,
        values=vals,
        origin=xmin,
        span=span,
        out_svg=out_svg,
        width_in=width_in,
        height_in=height_in,
    )


def export_consensus_zoom_manual_fixed_size(
    csv_path: str,
    res_start: int,
    res_end: int,
    out_svg: str,
    zoom_width_in: float = 1.6089,
    zoom_height_in: float = 0.0392,
    resi_col: str = "resi",
    val_col: str = "pctD_diff",
) -> None:
    """Export an exact-size Illustrator-ready zoom strip from the notebook."""
    if res_end < res_start:
        raise ValueError("res_end must be greater than or equal to res_start")

    df = _load_consensus(csv_path, resi_col, val_col)
    zoom_df = df[(df[resi_col] >= res_start) & (df[resi_col] <= res_end)].copy()
    if zoom_df.empty:
        raise ValueError(f"No rows found in range {res_start}-{res_end} using '{resi_col}'.")

    _write_manual_svg(
        residues=zoom_df[resi_col].to_numpy(dtype=int),
        values=zoom_df[val_col].to_numpy(),
        origin=int(res_start),
        span=int(res_end - res_start + 1),
        out_svg=out_svg,
        width_in=zoom_width_in,
        height_in=zoom_height_in,
    )


def export_svg_strip(
    csv_path: str,
    output: str,
    *,
    value_col: str = "pctD_diff",
    resi_col: str = "resi",
    resi_min: int | None = None,
    resi_max: int | None = None,
    width_in: float | None = None,
    height_in: float | None = None,
) -> None:
    """Backward-compatible dispatcher using only the Illustrator exporters.

    With no residue range, this is the notebook's full-strip exporter. With a
    range, it is the notebook's zoom exporter. Width/height default to the
    corresponding notebook values.
    """
    if resi_min is None and resi_max is None:
        export_consensus_full_manual_fixed_size(
            csv_path,
            output,
            width_in=7.4558 if width_in is None else width_in,
            height_in=0.0751 if height_in is None else height_in,
            resi_col=resi_col,
            val_col=value_col,
        )
        return

    if resi_min is None or resi_max is None:
        raise ValueError("Both resi_min and resi_max are required for a zoom strip")
    export_consensus_zoom_manual_fixed_size(
        csv_path,
        res_start=resi_min,
        res_end=resi_max,
        out_svg=output,
        zoom_width_in=1.6089 if width_in is None else width_in,
        zoom_height_in=0.0392 if height_in is None else height_in,
        resi_col=resi_col,
        val_col=value_col,
    )
