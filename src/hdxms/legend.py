"""Generate standalone HDX-MS color legends."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from hdxms.colors import CHIMERAX_HEX


LEGEND_ROWS = (
    ("hdx_blue3", "< -30", "Strong protection"),
    ("hdx_blue2", "-30 to < -15", "Moderate protection"),
    ("hdx_blue1", "-15 to < -5", "Weak protection"),
    ("hdx_gray0", "-5 to 5", "No significant change"),
    ("hdx_red1", "> 5 to 15", "Weak deprotection"),
    ("hdx_red2", "> 15 to 30", "Moderate deprotection"),
    ("hdx_red3", "> 30", "Strong deprotection"),
    ("hdx_unmapped", "NA", "No experimental coverage"),
)


def write_hdx_legend(
    *,
    svg_path: str | None = None,
    png_path: str | None = None,
    title: str = "HDX-MS Difference",
    value_label: str = "Delta deuteration (%)",
    dpi: int = 300,
) -> None:
    """Write a transparent standalone legend as SVG and/or PNG."""

    if not svg_path and not png_path:
        return

    fig, ax = plt.subplots(figsize=(4.8, 3.45))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, len(LEGEND_ROWS) + 1.35)
    ax.axis("off")

    ax.text(
        0.02,
        len(LEGEND_ROWS) + 1.03,
        title,
        fontsize=12,
        fontweight="bold",
        va="center",
    )
    ax.text(
        0.02,
        len(LEGEND_ROWS) + 0.48,
        value_label,
        fontsize=9,
        va="center",
    )

    top = len(LEGEND_ROWS) - 0.05
    for index, (color_name, interval, label) in enumerate(LEGEND_ROWS):
        y = top - index
        color = CHIMERAX_HEX[color_name]

        ax.add_patch(
            Rectangle(
                (0.02, y - 0.28),
                0.115,
                0.56,
                facecolor=color,
                edgecolor="#303030",
                linewidth=0.6,
            )
        )
        ax.text(
            0.17,
            y,
            interval,
            fontsize=8.5,
            family="DejaVu Sans Mono",
            va="center",
        )
        ax.text(
            0.48,
            y,
            label,
            fontsize=8.5,
            va="center",
        )

    fig.tight_layout(pad=0.35)

    for output_path, output_format in (
        (svg_path, "svg"),
        (png_path, "png"),
    ):
        if not output_path:
            continue

        output = Path(output_path).expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            output,
            format=output_format,
            dpi=dpi,
            bbox_inches="tight",
            transparent=True,
        )

    plt.close(fig)
