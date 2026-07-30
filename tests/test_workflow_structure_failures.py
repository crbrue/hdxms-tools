from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd
import pytest
import yaml

import hdxms.workflow as workflow


def _row(start, end, sequence, charge, pct):
    return {
        "Start": start,
        "End": end,
        "Sequence": sequence,
        "Charge": charge,
        "# Deut": pct / 25.0,
        "Deut %": pct,
        "Confidence": "High",
    }


def _manifest_project(tmp_path: Path, *, structure: bool = True, empirical: bool = False):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    datasets = {
        "fd.csv": [_row(1, 5, "AAAAA", 2, 80.0)],
        "a.csv": [_row(1, 5, "AAAAA", 2, 40.0)],
        "b.csv": [_row(1, 5, "AAAAA", 2, 20.0)],
    }
    for name, rows in datasets.items():
        pd.DataFrame(rows).to_csv(inputs / name, index=False)
    pdb = tmp_path / "protein.pdb"
    pdb.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n")
    comparison = {"name": "A_vs_B", "dataset_a": "A", "dataset_b": "B"}
    if structure:
        comparison["structure"] = "protein.pdb"
    manifest = {
        "base_dir": ".",
        "output_dir": "output",
        "fully_deuterated": "inputs/fd.csv",
        "sample_d2o_fraction": 0.75,
        "fd_d2o_fraction": 1.0,
        "apply_empirical_correction": empirical,
        "dpi": 50,
        "quality_control": {"enabled": False},
        "datasets": {
            "A": {"files": ["inputs/a.csv"]},
            "B": {"files": ["inputs/b.csv"]},
        },
        "comparisons": [comparison],
    }
    path = tmp_path / "manifest.yaml"
    path.write_text(yaml.safe_dump(manifest, sort_keys=False))
    return path


def _report_messages(output: Path):
    report = yaml.safe_load((output / "analysis_report/analysis_report.yaml").read_text())
    return [str(item.get("message", "")) for item in report.get("warnings", [])]


@pytest.mark.parametrize(
    "error, expected_fragment",
    [
        (FileNotFoundError("pymol"), "PyMOL executable not found"),
        (subprocess.CalledProcessError(9, ["pymol"]), "exit code 9"),
    ],
)
def test_pymol_failure_is_nonfatal_and_reported(tmp_path: Path, monkeypatch, error, expected_fragment):
    manifest = _manifest_project(tmp_path, structure=True)

    def fail(*args, **kwargs):
        raise error

    monkeypatch.setattr(workflow.subprocess, "run", fail)
    output = workflow.run_manifest(manifest)

    assert (output / "comparisons/A_vs_B/no_be/tables/pctD_summary.csv").is_file()
    assert (output / "comparisons/A_vs_B/be/tables/pctD_summary.csv").is_file()
    messages = _report_messages(output)
    assert sum(expected_fragment in message for message in messages) == 2


def test_manifest_without_structure_never_calls_pymol(tmp_path: Path, monkeypatch):
    manifest = _manifest_project(tmp_path, structure=False)

    def unexpected(*args, **kwargs):
        pytest.fail("PyMOL should not run without a structure entry")

    monkeypatch.setattr(workflow.subprocess, "run", unexpected)
    output = workflow.run_manifest(manifest)

    assert not (output / "comparisons/A_vs_B/no_be/structure").exists()
    assert not (output / "comparisons/A_vs_B/be/structure").exists()


def test_missing_structure_file_fails_before_pymol(tmp_path: Path, monkeypatch):
    manifest = _manifest_project(tmp_path, structure=True)
    data = yaml.safe_load(manifest.read_text())
    data["comparisons"][0]["structure"] = "does_not_exist.pdb"
    manifest.write_text(yaml.safe_dump(data, sort_keys=False))

    monkeypatch.setattr(workflow.subprocess, "run", lambda *a, **k: pytest.fail("PyMOL should not be called"))
    with pytest.raises(FileNotFoundError, match="does_not_exist.pdb"):
        workflow.run_manifest(manifest)


def test_empirical_mode_is_recorded_and_raw_branch_is_unchanged(tmp_path: Path, monkeypatch):
    manifest = _manifest_project(tmp_path, structure=False, empirical=True)
    monkeypatch.setattr(workflow.subprocess, "run", lambda *a, **k: None)
    output = workflow.run_manifest(manifest)

    no_be = pd.read_csv(output / "comparisons/A_vs_B/no_be/tables/pctD_summary.csv")
    be = pd.read_csv(output / "comparisons/A_vs_B/be/tables/pctD_summary.csv")
    assert no_be.loc[0, "pctD_diff"] == pytest.approx(20.0)
    assert set(no_be["normalization"]) == {"none"}
    assert set(be["normalization"]) == {"empirical"}
    metadata = yaml.safe_load((output / "run_metadata.yaml").read_text())
    assert metadata["sample_d2o_fraction"] == pytest.approx(0.75)
    assert metadata["fd_d2o_fraction"] == pytest.approx(1.0)
