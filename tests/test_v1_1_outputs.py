from pathlib import Path

import pandas as pd
import yaml

from hdxms.quality import build_back_exchange_table, build_peptide_quality_table
from hdxms.workflow import run_manifest


def _frame(values, pct_values=None):
    pct_values = pct_values or values
    return pd.DataFrame({
        "Start": [1, 10],
        "End": [5, 15],
        "Sequence": ["AAAAA", "AAPAAA"],
        "Charge": [2, 2],
        "# Deut": values,
        "Deut %": pct_values,
        "Confidence": ["High", "High"],
    })


def test_back_exchange_and_red_subset_logic():
    fd = _frame([2.0, 2.0], [50.0, 50.0])
    audit = build_back_exchange_table(fd, 0.75)
    raw = [_frame([1.0, 1.0]), _frame([1.1, 3.0]), _frame([0.9, 0.0])]
    corrected = [_frame([2.0, 2.0]), _frame([2.2, 8.0]), _frame([1.8, 0.0])]
    qc = build_peptide_quality_table("test", raw, corrected, audit, {
        "corrected_occupancy": {"warning": 1.0, "severe": 1.1},
        "replicate_sd_deuterons": {"warning": 0.5, "severe": 1.0},
    })
    assert {"QC Score", "QC Category", "QC Flags", "Recommendation"}.issubset(qc.columns)
    assert "FD_OVERCORRECTION_GT_110" in qc.loc[1, "QC Flags"]
    assert qc.loc[1, "QC Category"] == "Poor"


def test_manifest_creates_dual_si_and_report(tmp_path: Path):
    exports = tmp_path / "exports"
    exports.mkdir()
    fd = _frame([2.0, 2.0], [50.0, 50.0])
    fd.to_csv(exports / "fd.csv", index=False)
    for i, vals in enumerate(([1.0, 1.0], [1.1, 1.2], [0.9, 1.1]), 1):
        _frame(list(vals), [25.0, 25.0]).to_csv(exports / f"apo_rep{i}.csv", index=False)

    manifest = {
        "base_dir": ".",
        "output_dir": "out",
        "sample_d2o_fraction": 0.75,
        "fd_d2o_fraction": 0.75,
        "fully_deuterated": "exports/fd.csv",
        "datasets": {"apo": {"files": [f"exports/apo_rep{i}.csv" for i in (1, 2, 3)]}},
        "comparisons": [],
        "quality_control": {"enabled": True},
    }
    manifest_path = tmp_path / "manifest.yaml"
    manifest_path.write_text(yaml.safe_dump(manifest))
    out = run_manifest(manifest_path)

    no_be = pd.read_csv(out / "datasets/apo/no_be/SI_summary_table.csv")
    be = pd.read_csv(out / "datasets/apo/be/SI_summary_table.csv")
    assert list(no_be.columns) == list(be.columns)
    assert (out / "back_exchange/peptide_back_exchange.csv").exists()
    assert (out / "analysis_report/peptide_quality/peptide_quality_scores.csv").exists()
    assert (out / "analysis_report/peptide_quality/flagged_peptides.csv").exists()
    report = yaml.safe_load((out / "analysis_report/analysis_report.yaml").read_text())
    assert report["schema_version"] == "1.1"
    assert report["normalization"]["si_table_schema_preserved_between_modes"] is True
