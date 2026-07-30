"""MADP HDX Palette v1.0.

Canonical color definitions and binning helpers for MADP HDX-MS figures.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Optional


PALETTE_NAME = "MADP HDX Palette"
PALETTE_VERSION = "1.0"


@dataclass(frozen=True)
class HDXColor:
    id: str
    label: str
    hex: str
    rgb_255: tuple[int, int, int]
    rgb_fraction: tuple[float, float, float]


HDX_COLORS: dict[str, HDXColor] = {
    "hdx_blue3": HDXColor(
        "hdx_blue3", "Strong protection", "#3366E6",
        (51, 102, 230), (0.2, 0.4, 0.9),
    ),
    "hdx_blue2": HDXColor(
        "hdx_blue2", "Moderate protection", "#6699FF",
        (102, 153, 255), (0.4, 0.6, 1.0),
    ),
    "hdx_blue1": HDXColor(
        "hdx_blue1", "Weak protection", "#A6BFFF",
        (166, 191, 255), (166 / 255, 191 / 255, 1.0),
    ),
    "hdx_gray0": HDXColor(
        "hdx_gray0", "Neutral", "#595959",
        (89, 89, 89), (89 / 255, 89 / 255, 89 / 255),
    ),
    "hdx_red1": HDXColor(
        "hdx_red1", "Weak deprotection", "#FFB3B3",
        (255, 179, 179), (1.0, 179 / 255, 179 / 255),
    ),
    "hdx_red2": HDXColor(
        "hdx_red2", "Moderate deprotection", "#F27373",
        (242, 115, 115), (242 / 255, 115 / 255, 115 / 255),
    ),
    "hdx_red3": HDXColor(
        "hdx_red3", "Strong deprotection", "#D93333",
        (217, 51, 51), (217 / 255, 51 / 255, 51 / 255),
    ),
    "hdx_unmapped": HDXColor(
        "hdx_unmapped", "No coverage", "#A8A8A8",
        (168, 168, 168), (168 / 255, 168 / 255, 168 / 255),
    ),
}


def hdx_bin(value: float | int | None) -> str:
    """Return the canonical discrete color ID for an HDX value.

    Bin boundaries match the PyMOL/ChimeraX structure-coloring workflow:
      value < -30      -> hdx_blue3
      -30 <= value < -15 -> hdx_blue2
      -15 <= value <= -5 -> hdx_blue1
      -5 < value <= 5    -> hdx_gray0
      5 < value <= 15    -> hdx_red1
      15 < value <= 30   -> hdx_red2
      value > 30         -> hdx_red3
      missing/non-finite -> hdx_unmapped
    """
    if value is None:
        return "hdx_unmapped"

    numeric = float(value)
    if not isfinite(numeric):
        return "hdx_unmapped"
    if numeric < -30:
        return "hdx_blue3"
    if numeric < -15:
        return "hdx_blue2"
    if numeric <= -5:
        return "hdx_blue1"
    if numeric <= 5:
        return "hdx_gray0"
    if numeric <= 15:
        return "hdx_red1"
    if numeric <= 30:
        return "hdx_red2"
    return "hdx_red3"


def hex_for_value(value: float | int | None) -> str:
    return HDX_COLORS[hdx_bin(value)].hex


def rgb_fraction_for_value(
    value: float | int | None,
) -> tuple[float, float, float]:
    return HDX_COLORS[hdx_bin(value)].rgb_fraction


def pymol_color_definitions() -> list[tuple[str, tuple[float, float, float]]]:
    """Return entries suitable for cmd.set_color(name, rgb_fraction)."""
    return [(name, color.rgb_fraction) for name, color in HDX_COLORS.items()]
