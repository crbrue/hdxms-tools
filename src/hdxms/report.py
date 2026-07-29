from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


def write_analysis_report(
    output: str | Path,
    *,
    manifest: str,
    normalization: dict[str, Any],
    datasets: dict[str, Any],
    comparisons: list[dict[str, Any]],
    input_inventory_csv: str,
    peptide_attrition_csv: str,
    output_inventory: list[str],
    warnings: list[dict[str, Any]],
) -> Path:
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    report = {
        "schema_version": "1.1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": manifest,
        "normalization": normalization,
        "datasets": datasets,
        "comparisons": comparisons,
        "provenance": {
            "raw_input_inventory": input_inventory_csv,
            "peptide_attrition": peptide_attrition_csv,
        },
        "warnings": warnings,
        "output_inventory": sorted(output_inventory),
    }
    output.write_text(yaml.safe_dump(report, sort_keys=False))
    return output
