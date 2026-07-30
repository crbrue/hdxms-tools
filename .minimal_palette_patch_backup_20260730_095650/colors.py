"""Canonical MADP HDX color palette and binning rules."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable


PALETTE_NAME = "MADP HDX Palette"
PALETTE_VERSION = "1.0"


@dataclass(frozen=True)
class HDXColor:
    """One named color in the canonical MADP HDX palette."""

    name: str
    label: str
    hex: str
    rgb_255: tuple[int, int, int]

    @property
    def rgb_fraction(self) -> tuple[float, float, float]:
        """Return RGB channels normalized to the PyMOL 0–1 range."""
        return tuple(channel / 255.0 for channel in self.rgb_255)


COLORS: dict[str, HDXColor] = {
    "hdx_blue3": HDXColor(
        "hdx_blue3",
        "Strong protection",
        "#3366E6",
        (51, 102, 230),
    ),
    "hdx_blue2": HDXColor(
        "hdx_blue2",
        "Moderate protection",
        "#6699FF",
        (102, 153, 255),
    ),
    "hdx_blue1": HDXColor(
        "hdx_blue1",
        "Weak protection",
        "#A6BFFF",
        (166, 191, 255),
    ),
    "hdx_gray0": HDXColor(
        "hdx_gray0",
        "Neutral",
        "#595959",
        (89, 89, 89),
    ),
    "hdx_red1": HDXColor(
        "hdx_red1",
        "Weak deprotection",
        "#FFB3B3",
        (255, 179, 179),
    ),
    "hdx_red2": HDXColor(
        "hdx_red2",
        "Moderate deprotection",
        "#F27373",
        (242, 115, 115),
    ),
    "hdx_red3": HDXColor(
        "hdx_red3",
        "Strong deprotection",
        "#D93333",
        (217, 51, 51),
    ),
    "hdx_unmapped": HDXColor(
        "hdx_unmapped",
        "No experimental coverage",
        "#A8A8A8",
        (168, 168, 168),
    ),
}


# Preserve the original hdxms.colors API used by diffplots.py and potentially
# by downstream user scripts. Each entry is:
#
#     (PyMOL color name, hexadecimal color, membership predicate)
#
# Boundary behavior is unchanged from the pre-refactor implementation.
HDX_BINS: tuple[
    tuple[str, str, Callable[[float], bool]],
    ...,
] = (
    (
        "hdx_blue3",
        COLORS["hdx_blue3"].hex,
        lambda value: value < -30,
    ),
    (
        "hdx_blue2",
        COLORS["hdx_blue2"].hex,
        lambda value: -30 <= value < -15,
    ),
    (
        "hdx_blue1",
        COLORS["hdx_blue1"].hex,
        lambda value: -15 <= value < -5,
    ),
    (
        "hdx_gray0",
        COLORS["hdx_gray0"].hex,
        lambda value: -5 <= value <= 5,
    ),
    (
        "hdx_red1",
        COLORS["hdx_red1"].hex,
        lambda value: 5 < value <= 15,
    ),
    (
        "hdx_red2",
        COLORS["hdx_red2"].hex,
        lambda value: 15 < value <= 30,
    ),
    (
        "hdx_red3",
        COLORS["hdx_red3"].hex,
        lambda value: value > 30,
    ),
)


# PyMOL consumes normalized RGB values.
PYMOL_RGB: dict[str, list[float]] = {
    name: list(color.rgb_fraction)
    for name, color in COLORS.items()
}


# ChimeraX and standard plotting code consume hexadecimal values.
CHIMERAX_HEX: dict[str, str] = {
    name: color.hex
    for name, color in COLORS.items()
}


# General backwards-compatible alias.
HEX_COLORS: dict[str, str] = CHIMERAX_HEX.copy()


def bin_name(value: float | int | None) -> str | None:
    """Return the canonical discrete HDX color-bin name.

    Missing, NaN, and infinite values return ``None``.
    """

    if value is None:
        return None

    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None

    if not math.isfinite(numeric):
        return None

    for name, _, predicate in HDX_BINS:
        if predicate(numeric):
            return name

    return None


def hex_color(
    value: float | int | None,
    missing: str = "#ffffff",
) -> str:
    """Return the canonical hexadecimal color for an HDX value.

    The default missing color remains white for compatibility with the
    existing difference-map plotting functions. Callers may explicitly pass
    ``COLORS["hdx_unmapped"].hex`` when no-coverage residues should appear
    gray.
    """

    name = bin_name(value)

    if name is None:
        return missing

    return COLORS[name].hex


def rgb_fraction(
    value: float | int | None,
) -> tuple[float, float, float]:
    """Return a normalized RGB tuple for an HDX value.

    Missing and non-finite values use the no-coverage color.
    """

    name = bin_name(value) or "hdx_unmapped"
    return COLORS[name].rgb_fraction
