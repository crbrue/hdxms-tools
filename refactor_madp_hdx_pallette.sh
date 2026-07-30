#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${1:-$PWD}"
PKG_DIR="${REPO_ROOT}/src/hdxms"

for required in \
    "${PKG_DIR}/colors.py" \
    "${PKG_DIR}/pymol_script.py" \
    "${PKG_DIR}/cli.py"
do
    if [[ ! -f "${required}" ]]; then
        echo "ERROR: missing ${required}" >&2
        echo "Run this from the hdx-ms-tools repository root, or pass the repo path." >&2
        exit 1
    fi
done

STAMP="$(date +%Y%m%d_%H%M%S)"
BACKUP_DIR="${REPO_ROOT}/.palette_refactor_backup_${STAMP}"
mkdir -p "${BACKUP_DIR}"

cp -v "${PKG_DIR}/colors.py" "${BACKUP_DIR}/colors.py"
cp -v "${PKG_DIR}/pymol_script.py" "${BACKUP_DIR}/pymol_script.py"
cp -v "${PKG_DIR}/cli.py" "${BACKUP_DIR}/cli.py"
[[ -f "${PKG_DIR}/workflow.py" ]] && cp -v "${PKG_DIR}/workflow.py" "${BACKUP_DIR}/workflow.py"

cat > "${PKG_DIR}/colors.py" <<'PY'
"""Canonical MADP HDX color palette and binning rules."""

from __future__ import annotations

import math
from dataclasses import dataclass


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
        return tuple(channel / 255.0 for channel in self.rgb_255)


COLORS: dict[str, HDXColor] = {
    "hdx_blue3": HDXColor(
        "hdx_blue3", "Strong protection", "#3366E6", (51, 102, 230)
    ),
    "hdx_blue2": HDXColor(
        "hdx_blue2", "Moderate protection", "#6699FF", (102, 153, 255)
    ),
    "hdx_blue1": HDXColor(
        "hdx_blue1", "Weak protection", "#A6BFFF", (166, 191, 255)
    ),
    "hdx_gray0": HDXColor(
        "hdx_gray0", "Neutral", "#595959", (89, 89, 89)
    ),
    "hdx_red1": HDXColor(
        "hdx_red1", "Weak deprotection", "#FFB3B3", (255, 179, 179)
    ),
    "hdx_red2": HDXColor(
        "hdx_red2", "Moderate deprotection", "#F27373", (242, 115, 115)
    ),
    "hdx_red3": HDXColor(
        "hdx_red3", "Strong deprotection", "#D93333", (217, 51, 51)
    ),
    "hdx_unmapped": HDXColor(
        "hdx_unmapped", "No coverage", "#A8A8A8", (168, 168, 168)
    ),
}

# Backwards-compatible exports used elsewhere in hdx-ms-tools.
PYMOL_RGB: dict[str, list[float]] = {
    name: list(color.rgb_fraction)
    for name, color in COLORS.items()
}

CHIMERAX_HEX: dict[str, str] = {
    name: color.hex
    for name, color in COLORS.items()
}

HEX_COLORS: dict[str, str] = CHIMERAX_HEX.copy()


def bin_name(value: float | int | None) -> str | None:
    """Return the canonical discrete HDX bin for a numerical value.

    Boundaries:
        value < -30       strong protection
        -30 <= value < -15
        -15 <= value <= -5
        -5 < value <= 5   neutral
        5 < value <= 15
        15 < value <= 30
        value > 30        strong deprotection

    Missing and non-finite values return ``None``.
    """
    if value is None:
        return None

    numeric = float(value)
    if not math.isfinite(numeric):
        return None
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


def hex_color(value: float | int | None, missing: str | None = None) -> str:
    """Return the canonical hex color for a value."""
    name = bin_name(value)
    if name is None:
        return missing if missing is not None else COLORS["hdx_unmapped"].hex
    return COLORS[name].hex


def rgb_fraction(value: float | int | None) -> tuple[float, float, float]:
    """Return a normalized RGB tuple for a value."""
    name = bin_name(value) or "hdx_unmapped"
    return COLORS[name].rgb_fraction
PY

cat > "${PKG_DIR}/pymol_script.py" <<'PY'
from __future__ import annotations

import argparse
import math
import os
from collections import defaultdict
from pathlib import Path

from pymol import cmd, util

from hdxms.colors import CHIMERAX_HEX, PYMOL_RGB, bin_name
from hdxms.io import apply_offsets, load_hdx_csv


def _ranges(values):
    """Collapse residue numbers into compact contiguous ranges."""
    values = sorted(set(int(value) for value in values))
    if not values:
        return []

    output = []
    start = end = values[0]

    for value in values[1:]:
        if value == end + 1:
            end = value
        else:
            output.append(f"{start}-{end}" if start != end else str(start))
            start = end = value

    output.append(f"{start}-{end}" if start != end else str(start))
    return output


def _selection(obj: str, chain: str, residues) -> str:
    residue_expression = "+".join(_ranges(residues))
    chain_clause = f" and chain {chain}" if chain else ""
    return (
        f"{obj}{chain_clause} and polymer.protein "
        f"and resi {residue_expression}"
    )


def _write_chimerax_script(
    output_path: str,
    structure_path: str,
    grouped,
) -> None:
    """Write a ChimeraX script using the exact canonical discrete palette."""
    output = Path(output_path).expanduser().resolve()
    structure = Path(structure_path).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "# MADP HDX Palette v1.0",
        "# Generated by hdxms color",
        f'open "{structure}"',
        "hide",
        "show protein cartoons",
        f"color protein {CHIMERAX_HEX['hdx_unmapped']} target c",
        "show nucleic atoms",
        "style nucleic stick",
        "color nucleic byelement",
        "show ligand atoms",
        "style ligand stick",
        "color ligand byelement",
        "",
        "# Discrete HDX residue colors",
    ]

    color_order = (
        "hdx_blue3",
        "hdx_blue2",
        "hdx_blue1",
        "hdx_gray0",
        "hdx_red1",
        "hdx_red2",
        "hdx_red3",
    )

    for color_name in color_order:
        for chain, residues in sorted(grouped.get(color_name, {}).items()):
            chain_prefix = f"/{chain}" if chain else ""
            for residue_range in _ranges(residues):
                lines.append(
                    f"color {chain_prefix}:{residue_range} "
                    f"{CHIMERAX_HEX[color_name]} target c"
                )

    lines.extend(
        [
            "",
            "set bgColor white",
            "lighting soft",
            "view protein",
        ]
    )

    output.write_text("\n".join(lines) + "\n")
    print(f"Wrote ChimeraX script: {output}")


def run(args):
    obj = os.path.splitext(os.path.basename(args.structure))[0]
    cmd.load(args.structure, obj)

    for name, rgb in PYMOL_RGB.items():
        cmd.set_color(name, rgb)

    df = apply_offsets(
        load_hdx_csv(args.csv, args.value_col),
        args.offset,
        args.offset_per_chain,
    )

    allowed = None
    if args.chain_filter:
        allowed = {
            chain.strip()
            for chain in args.chain_filter.split(",")
            if chain.strip()
        }
        df = df[(df.chain == "") | df.chain.isin(allowed)]

    sums = defaultdict(float)
    counts = defaultdict(int)

    for row in df.itertuples(index=False):
        value = float(row.value)
        if not math.isfinite(value):
            continue

        for residue in range(int(row.start), int(row.end) + 1):
            key = (row.chain, residue)
            sums[key] += value
            counts[key] += 1

    model = cmd.get_model(f"{obj} and polymer.protein")
    present = defaultdict(set)

    for atom in model.atom:
        try:
            present[atom.chain].add(int(atom.resi))
        except ValueError:
            continue

    mapped = {}
    for chain, residues in present.items():
        if allowed is not None and chain not in allowed:
            continue

        for residue in residues:
            chain_key = (chain, residue)
            blank_key = ("", residue)
            key = chain_key if chain_key in sums else blank_key

            if counts.get(key, 0) >= args.mincov:
                mapped[(chain, residue)] = sums[key] / counts[key]

    cmd.hide("everything", "all")
    cmd.show("cartoon", f"{obj} and polymer.protein")
    cmd.color("hdx_unmapped", f"{obj} and polymer.protein")

    grouped = defaultdict(lambda: defaultdict(list))
    for (chain, residue), value in mapped.items():
        color_name = bin_name(value)
        if color_name is not None:
            grouped[color_name][chain].append(residue)

    for color_name, chain_map in grouped.items():
        for chain, residues in chain_map.items():
            cmd.color(color_name, _selection(obj, chain, residues))

    if not args.protein_only:
        cmd.show("sticks", f"{obj} and polymer.nucleic")
        util.cbaw(f"{obj} and polymer.nucleic")
        cmd.show("sticks", f"{obj} and organic")
        util.cbag(f"{obj} and organic")

    cmd.set("cartoon_transparency", args.cartoon_trans)
    cmd.bg_color("white")
    cmd.orient(f"{obj} and polymer.protein")

    if args.out_bfactor_pdb:
        # Preserve original B factors for unmapped atoms and update only residues
        # that carry finite HDX values.
        for (chain, residue), value in mapped.items():
            chain_clause = f" and chain {chain}" if chain else ""
            cmd.alter(
                f"{obj}{chain_clause} and polymer.protein and resi {residue}",
                f"b={float(value)}",
            )
        cmd.save(args.out_bfactor_pdb, obj)

    if args.out_chimerax:
        chimerax_structure = args.out_bfactor_pdb or args.structure
        _write_chimerax_script(
            args.out_chimerax,
            chimerax_structure,
            grouped,
        )

    if args.session:
        cmd.save(args.session)

    if args.out_png:
        cmd.viewport(args.width, args.height)
        cmd.ray()
        cmd.png(args.out_png, dpi=args.dpi)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("structure")
    parser.add_argument("csv")
    parser.add_argument("--value-col", default="pctD_diff")
    parser.add_argument("--chain-filter")
    parser.add_argument("--offset", type=int, default=0)
    parser.add_argument("--offset-per-chain")
    parser.add_argument("--mincov", type=int, default=1)
    parser.add_argument("--cartoon-trans", type=float, default=0.0)
    parser.add_argument("--protein-only", action="store_true")
    parser.add_argument("--out-png")
    parser.add_argument("--session")
    parser.add_argument("--out-bfactor-pdb")
    parser.add_argument("--out-chimerax")
    parser.add_argument("--width", type=int, default=1800)
    parser.add_argument("--height", type=int, default=1400)
    parser.add_argument("--dpi", type=int, default=300)
    run(parser.parse_args(argv))


if __name__ == "__main__":
    main()
PY

python - "${PKG_DIR}/cli.py" "${PKG_DIR}/workflow.py" <<'PY'
from pathlib import Path
import re
import sys

cli_path = Path(sys.argv[1])
workflow_path = Path(sys.argv[2])

text = cli_path.read_text()

# Add the public CLI option if absent.
if "--out-chimerax" not in text:
    target = "p.add_argument('--out-bfactor-pdb')"
    if target not in text:
        raise SystemExit("ERROR: could not find --out-bfactor-pdb in cli.py")
    text = text.replace(
        target,
        target + "; p.add_argument('--out-chimerax')",
        1,
    )

# Ensure the option is forwarded to pymol_script.py.
if "('--out-chimerax',args.out_chimerax)" not in text:
    pattern = re.compile(
        r"\["
        r"\('--chain-filter',args\.chain_filter\),"
        r"\('--offset-per-chain',args\.offset_per_chain\),"
        r"\('--out-png',args\.out_png\),"
        r"\('--session',args\.session\),"
        r"\('--out-bfactor-pdb',args\.out_bfactor_pdb\)"
        r"\]"
    )
    replacement = (
        "[('--chain-filter',args.chain_filter),"
        "('--offset-per-chain',args.offset_per_chain),"
        "('--out-png',args.out_png),"
        "('--session',args.session),"
        "('--out-bfactor-pdb',args.out_bfactor_pdb),"
        "('--out-chimerax',args.out_chimerax)]"
    )
    text, count = pattern.subn(replacement, text, count=1)
    if count != 1:
        raise SystemExit(
            "ERROR: could not find the color output-forwarding list in cli.py"
        )

# PyMOL must explicitly run the Python file.
text = text.replace(
    "cmd=[args.pymol,'-cq',str(module),'--'",
    "cmd=[args.pymol,'-cq','-r',str(module),'--'",
)
text = text.replace(
    'cmd=[args.pymol, "-cq", str(module), "--"',
    'cmd=[args.pymol, "-cq", "-r", str(module), "--"',
)

if "'-r',str(module)" not in text and '"-r", str(module)' not in text:
    raise SystemExit("ERROR: could not confirm the cli.py PyMOL -r fix")

compile(text, str(cli_path), "exec")
cli_path.write_text(text)

# Keep manifest-driven rendering consistent where workflow.py exists.
if workflow_path.exists():
    workflow = workflow_path.read_text()
    workflow = workflow.replace(
        '"-cq",\n                    str(pymol_module),',
        '"-cq",\n                    "-r",\n                    str(pymol_module),',
    )
    workflow = workflow.replace(
        "'-cq',\n                    str(pymol_module),",
        "'-cq',\n                    '-r',\n                    str(pymol_module),",
    )

    # Add a ChimeraX product next to the other structure products.
    if '"--out-chimerax"' not in workflow and "'--out-chimerax'" not in workflow:
        marker = (
            '"--out-bfactor-pdb", '
            'str(structure_dir / "consensus_structure_bfactor.pdb"),'
        )
        if marker in workflow:
            workflow = workflow.replace(
                marker,
                marker
                + '\n                    "--out-chimerax", '
                'str(structure_dir / "consensus_structure.cxc"),',
                1,
            )

    compile(workflow, str(workflow_path), "exec")
    workflow_path.write_text(workflow)
PY

python -m compileall -q "${PKG_DIR}"

echo
echo "Refactor installed successfully."
echo "Backup directory:"
echo "  ${BACKUP_DIR}"
echo
echo "Verify:"
echo "  hdxms color --help | grep -E 'out-chimerax|out-bfactor'"
echo
echo "Recommended test:"
echo "  python -m pytest -q tests/test_color_cli.py tests/test_pymol_script.py"