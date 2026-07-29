#!/usr/bin/env python3
"""Create copy/paste-ready CSV filename inventories for an HDX-MS manifest.

The primary output is deliberately inference-free: every CSV filename in the
raw-data directory is written to a Python list.  Optional manifest inference is
best-effort and never removes or silently assigns ambiguous files.
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

FD_RE = re.compile(r"(^|[_\-.])(fd|full[_\-.]?d|fully[_\-.]?deuterated)([_\-.]|$)", re.I)
REP_RE = re.compile(r"(?:[_\-.](?:rep|r)?\d+)$", re.I)


def csv_names(raw_dir: Path, recursive: bool = False) -> list[str]:
    pattern = "**/*.csv" if recursive else "*.csv"
    files = [p for p in raw_dir.glob(pattern) if p.is_file()]
    return sorted((p.relative_to(raw_dir).as_posix() for p in files), key=str.casefold)


def python_list_text(names: list[str], variable: str = "raw_files") -> str:
    lines = [f"{variable} = ["]
    lines.extend(f"    {name!r}," for name in names)
    lines.append("]")
    return "\n".join(lines) + "\n"


def normalized_group_stem(name: str) -> str:
    stem = Path(name).stem
    stem = FD_RE.sub("_", stem)
    stem = REP_RE.sub("", stem)
    return re.sub(r"[_\-.]+", "_", stem).strip("_")


def infer_groups(names: list[str]) -> tuple[dict[str, dict[str, object]], list[str]]:
    """Infer dataset groups conservatively.

    A group is emitted only when it has at least one non-FD CSV.  Exactly one FD
    candidate is accepted. Zero or multiple FD candidates remain explicit in
    the output so the user must resolve them.
    """
    buckets: dict[str, dict[str, list[str]]] = defaultdict(lambda: {"raw": [], "fd": []})
    for name in names:
        key = normalized_group_stem(name)
        if not key:
            key = "UNRESOLVED"
        bucket = buckets[key]
        (bucket["fd"] if FD_RE.search(Path(name).stem) else bucket["raw"]).append(name)

    groups: dict[str, dict[str, object]] = {}
    ambiguous: list[str] = []
    for key, bucket in sorted(buckets.items(), key=lambda item: item[0].casefold()):
        raw = sorted(bucket["raw"], key=str.casefold)
        fd = sorted(bucket["fd"], key=str.casefold)
        if not raw:
            ambiguous.extend(fd)
            continue
        entry: dict[str, object] = {"raw_files": raw}
        if len(fd) == 1:
            entry["fd_file"] = fd[0]
        elif fd:
            entry["fd_candidates"] = fd
            ambiguous.extend(fd)
        groups[key] = entry
    return groups, sorted(set(ambiguous), key=str.casefold)


def yaml_suggestion(names: list[str]) -> str:
    groups, ambiguous = infer_groups(names)
    lines = [
        "# BEST-EFFORT SUGGESTION — REVIEW BEFORE USING",
        "# The complete, authoritative filename list is in raw_csv_files.py.",
        "datasets:",
    ]
    for name, entry in groups.items():
        lines.append(f"  - name: {name}")
        lines.append("    raw_files:")
        for filename in entry["raw_files"]:  # type: ignore[index]
            lines.append(f"      - {filename}")
        if "fd_file" in entry:
            lines.append(f"    fd_file: {entry['fd_file']}")
        elif "fd_candidates" in entry:
            lines.append("    # Multiple FD candidates detected; choose one:")
            for filename in entry["fd_candidates"]:  # type: ignore[index]
                lines.append(f"    #   - {filename}")
        else:
            lines.append("    # fd_file: ADD_IF_APPLICABLE")
    if ambiguous:
        lines.append("")
        lines.append("# Files requiring manual assignment:")
        lines.extend(f"#   - {filename}" for filename in ambiguous)
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="List raw-data CSV files and optionally suggest manifest groupings."
    )
    parser.add_argument("raw_dir", type=Path, help="Directory containing raw CSV files")
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("raw_csv_files.py"),
        help="Python-list output file (default: raw_csv_files.py)",
    )
    parser.add_argument(
        "--variable", default="raw_files", help="Python variable name (default: raw_files)"
    )
    parser.add_argument("--recursive", action="store_true", help="Include CSV files in subdirectories")
    parser.add_argument(
        "--suggest-manifest", type=Path, metavar="FILE",
        help="Also write a best-effort YAML dataset skeleton",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    if not raw_dir.is_dir():
        raise SystemExit(f"ERROR: raw-data directory does not exist: {raw_dir}")

    names = csv_names(raw_dir, recursive=args.recursive)
    if not names:
        raise SystemExit(f"ERROR: no .csv files found in: {raw_dir}")

    args.output.write_text(python_list_text(names, args.variable), encoding="utf-8")
    print(f"Wrote {len(names)} CSV filenames to {args.output}")

    if args.suggest_manifest:
        args.suggest_manifest.write_text(yaml_suggestion(names), encoding="utf-8")
        print(f"Wrote best-effort manifest suggestion to {args.suggest_manifest}")
        print("Review inferred groupings before copying them into the manifest.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
