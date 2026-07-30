# MADP HDX Palette v1.0

Canonical discrete color palette for MADP HDX-MS consensus visualizations.

## Interpretation

| ID | Meaning | HDX interval | Hex | RGB |
|---|---|---:|---|---:|
| `hdx_blue3` | Strong protection | `< -30` | `#3366E6` | 51, 102, 230 |
| `hdx_blue2` | Moderate protection | `-30 to < -15` | `#6699FF` | 102, 153, 255 |
| `hdx_blue1` | Weak protection | `-15 to -5` | `#A6BFFF` | 166, 191, 255 |
| `hdx_gray0` | Neutral | `-5 to +5` | `#595959` | 89, 89, 89 |
| `hdx_red1` | Weak deprotection | `> +5 to +15` | `#FFB3B3` | 255, 179, 179 |
| `hdx_red2` | Moderate deprotection | `> +15 to +30` | `#F27373` | 242, 115, 115 |
| `hdx_red3` | Strong deprotection | `> +30` | `#D93333` | 217, 51, 51 |
| `hdx_unmapped` | No experimental coverage | N/A | `#A8A8A8` | 168, 168, 168 |

Negative differences indicate protection and positive differences indicate
deprotection. Non-finite and missing values must be assigned
`hdx_unmapped`, not the neutral bin.

## Boundary rules

The implementation uses these exact comparisons:

```python
if value < -30:
    hdx_blue3
elif value < -15:
    hdx_blue2
elif value <= -5:
    hdx_blue1
elif value <= 5:
    hdx_gray0
elif value <= 15:
    hdx_red1
elif value <= 30:
    hdx_red2
else:
    hdx_red3
```

## Repository usage

The canonical source of truth is `hdx_palette.yaml`.

Python:

```python
from hdxms.styles.hdx_palette import HDX_COLORS, hdx_bin
color = HDX_COLORS[hdx_bin(value)]
```

PyMOL:

```python
from resources.styles.pymol_colors import MADP_HDX_COLORS
for name, rgb in MADP_HDX_COLORS.items():
    cmd.set_color(name, rgb)
```

ChimeraX:

```text
color /A:100-110 #3366E6 target c
color protein #A8A8A8 target c
```

## Versioning

Increment the palette version whenever a bin boundary, RGB/hex value, or
semantic meaning changes. Figure-generation metadata should record the palette
name and version.
