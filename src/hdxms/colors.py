from __future__ import annotations
import math

HDX_BINS = (
    ("hdx_blue3", "#282a73", lambda v: v < -30),
    ("hdx_blue2", "#4478bb", lambda v: -30 <= v < -15),
    ("hdx_blue1", "#bde6f5", lambda v: -15 <= v < -5),
    ("hdx_gray0", "#afb0b0", lambda v: -5 <= v <= 5),
    ("hdx_red1", "#fac9c9", lambda v: 5 < v <= 15),
    ("hdx_red2", "#ed2024", lambda v: 15 < v <= 30),
    ("hdx_red3", "#7f1416", lambda v: v > 30),
)

PYMOL_RGB = {
    "hdx_blue3": [0.157, 0.165, 0.451],
    "hdx_blue2": [0.267, 0.471, 0.733],
    "hdx_blue1": [0.741, 0.902, 0.961],
    "hdx_gray0": [0.686, 0.690, 0.690],
    "hdx_red1": [0.980, 0.788, 0.788],
    "hdx_red2": [0.929, 0.125, 0.141],
    "hdx_red3": [0.498, 0.078, 0.086],
    "hdx_unmapped": [0.75, 0.75, 0.75],
}

def bin_name(value: float) -> str | None:
    if value is None or math.isnan(float(value)):
        return None
    value = float(value)
    for name, _, pred in HDX_BINS:
        if pred(value):
            return name
    return None

def hex_color(value: float, missing: str = "#ffffff") -> str:
    name = bin_name(value)
    if name is None:
        return missing
    return next(hex_value for n, hex_value, _ in HDX_BINS if n == name)
