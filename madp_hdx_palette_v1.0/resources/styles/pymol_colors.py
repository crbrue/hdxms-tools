"""Register MADP HDX Palette v1.0 in PyMOL."""

from pymol import cmd

MADP_HDX_COLORS = {
    "hdx_blue3": [0.2, 0.4, 0.9],
    "hdx_blue2": [0.4, 0.6, 1.0],
    "hdx_blue1": [166 / 255, 191 / 255, 1.0],
    "hdx_gray0": [89 / 255, 89 / 255, 89 / 255],
    "hdx_red1": [1.0, 179 / 255, 179 / 255],
    "hdx_red2": [242 / 255, 115 / 255, 115 / 255],
    "hdx_red3": [217 / 255, 51 / 255, 51 / 255],
    "hdx_unmapped": [168 / 255, 168 / 255, 168 / 255],
}

for color_name, rgb in MADP_HDX_COLORS.items():
    cmd.set_color(color_name, rgb)
