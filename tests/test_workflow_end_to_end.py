from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd
import pytest
import yaml

import hdxms.workflow as workflow


def _write_hx_csv(path: Path, rows: list[dict]) -> None:
    """Write the smallest HX Examiner-like CSV accepted by run_manifest()."""
    pd.DataFrame(rows).to_csv(path, index=False)


def _row(
    start: int,
    end: int,
    sequence: str,
    charge: int,
    deut_percent: float,
    deut_count: float,
) -> dict:
    return {
        "Start": start,
        "End": end,
        "Sequence": sequence,
        "Charge": charge,
        "# Deut": deut_count,
        "Deut %": deut_percent,
        "Confidence": "High",
    }


def _fake_pymol_run(calls: list[list[str]], real_run):
    """Intercept PyMOL renders while delegating unrelated subprocess calls.

    ``platform.platform()`` may internally call ``subprocess.run(["uname", "-p"])``,
    so the workflow-level monkeypatch must not treat every subprocess as PyMOL.
    """

    def fake_run(cmd, check=False, **kwargs):
        command = [str(value) for value in cmd]

        if "--out-png" not in command:
            return real_run(cmd, check=check, **kwargs)

        calls.append(command)
        assert check is True

        # Mimic the files that pymol_script.py would create.
        for flag in ("--out-png", "--session", "--out-bfactor-pdb"):
            output = Path(command[command.index(flag) + 1])
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(f"mock output for {flag}\n")

        return subprocess.CompletedProcess(command, 0)

    return fake_run


def test_manifest_workflow_generates_dual_comparisons_and_structure_outputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """Run a complete synthetic analysis and verify its numerical/output contract."""
    inputs = tmp_path / "inputs"
    inputs.mkdir()

    # FD is 80% D for both peptide identities.
    fd_rows = [
        _row(1, 5, "AAAAA", 2, 80.0, 3.20),
        _row(6, 10, "VVVVV", 3, 80.0, 3.20),
    ]

    # Raw replicate means:
    # AAAAA: A=32, B=12 -> no_be difference = +20
    # VVVVV: A=62, B=52 -> no_be difference = +10
    # FD normalization by 80% scales those differences by 100/80:
    # AAAAA -> +25; VVVVV -> +12.5
    a1_rows = [
        _row(1, 5, "AAAAA", 2, 30.0, 1.20),
        _row(6, 10, "VVVVV", 3, 60.0, 2.40),
    ]
    a2_rows = [
        _row(1, 5, "AAAAA", 2, 34.0, 1.36),
        _row(6, 10, "VVVVV", 3, 64.0, 2.56),
    ]
    b1_rows = [
        _row(1, 5, "AAAAA", 2, 10.0, 0.40),
        _row(6, 10, "VVVVV", 3, 50.0, 2.00),
    ]
    b2_rows = [
        _row(1, 5, "AAAAA", 2, 14.0, 0.56),
        _row(6, 10, "VVVVV", 3, 54.0, 2.16),
    ]

    files = {
        "fd.csv": fd_rows,
        "a_rep1.csv": a1_rows,
        "a_rep2.csv": a2_rows,
        "b_rep1.csv": b1_rows,
        "b_rep2.csv": b2_rows,
    }
    for filename, rows in files.items():
        _write_hx_csv(inputs / filename, rows)

    structure = tmp_path / "protein.pdb"
    structure.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )

    manifest = {
        "base_dir": ".",
        "output_dir": "output",
        "fully_deuterated": "inputs/fd.csv",
        "sample_d2o_fraction": 0.80,
        "fd_d2o_fraction": 0.80,
        "apply_empirical_correction": False,
        "dpi": 72,
        "quality_control": {"enabled": False},
        "datasets": {
            "A": {
                "title": "Condition A",
                "files": ["inputs/a_rep1.csv", "inputs/a_rep2.csv"],
            },
            "B": {
                "title": "Condition B",
                "files": ["inputs/b_rep1.csv", "inputs/b_rep2.csv"],
            },
        },
        "comparisons": [
            {
                "name": "A_vs_B",
                "dataset_a": "A",
                "dataset_b": "B",
                "label": "Condition A - Condition B",
                "structure": "protein.pdb",
                "chain_filter": "A",
                "protein_only": True,
            }
        ],
    }
    manifest_path = tmp_path / "manifest.yaml"
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False))

    pymol_calls: list[list[str]] = []
    real_subprocess_run = subprocess.run
    monkeypatch.setattr(
        workflow.subprocess,
        "run",
        _fake_pymol_run(pymol_calls, real_subprocess_run),
    )

    output = workflow.run_manifest(manifest_path)
    assert output == tmp_path / "output"

    expected_files = [
        # Dataset-level SI outputs.
        "datasets/A/no_be/SI_summary_table.csv",
        "datasets/A/be/SI_summary_table.csv",
        "datasets/B/no_be/SI_summary_table.csv",
        "datasets/B/be/SI_summary_table.csv",
        # Dual comparison branches.
        "comparisons/A_vs_B/no_be/tables/pctD_summary.csv",
        "comparisons/A_vs_B/no_be/tables/consensus_residue.csv",
        "comparisons/A_vs_B/no_be/figures/pctD_diff_bars.png",
        "comparisons/A_vs_B/no_be/figures/peptide_diff_map_stepwise.png",
        "comparisons/A_vs_B/no_be/figures/peptide_diff_map_compact.png",
        "comparisons/A_vs_B/no_be/figures/consensus_domain_diffmap.png",
        "comparisons/A_vs_B/no_be/illustrator/consensus_full.svg",
        "comparisons/A_vs_B/be/tables/pctD_summary.csv",
        "comparisons/A_vs_B/be/tables/consensus_residue.csv",
        "comparisons/A_vs_B/be/figures/pctD_diff_bars.png",
        "comparisons/A_vs_B/be/figures/peptide_diff_map_stepwise.png",
        "comparisons/A_vs_B/be/figures/peptide_diff_map_compact.png",
        "comparisons/A_vs_B/be/figures/consensus_domain_diffmap.png",
        "comparisons/A_vs_B/be/illustrator/consensus_full.svg",
        # Mocked PyMOL products for both branches.
        "comparisons/A_vs_B/no_be/structure/consensus_structure.png",
        "comparisons/A_vs_B/no_be/structure/consensus_structure.pse",
        "comparisons/A_vs_B/no_be/structure/consensus_structure_bfactor.pdb",
        "comparisons/A_vs_B/be/structure/consensus_structure.png",
        "comparisons/A_vs_B/be/structure/consensus_structure.pse",
        "comparisons/A_vs_B/be/structure/consensus_structure_bfactor.pdb",
        # Provenance/reporting contract.
        "provenance/manifest_used.yaml",
        "provenance/raw_input_inventory.csv",
        "provenance/peptide_attrition_by_input.csv",
        "provenance/peptide_attrition_by_dataset.csv",
        "provenance/software_environment.yaml",
        "analysis_report/analysis_report.yaml",
        "run_metadata.yaml",
    ]
    missing = [relative for relative in expected_files if not (output / relative).is_file()]
    assert missing == []

    no_be = pd.read_csv(output / "comparisons/A_vs_B/no_be/tables/pctD_summary.csv")
    be = pd.read_csv(output / "comparisons/A_vs_B/be/tables/pctD_summary.csv")

    # Peptide identity and charge are preserved through both branches.
    identity_columns = ["Start", "End", "Sequence", "Charge"]
    pd.testing.assert_frame_equal(no_be[identity_columns], be[identity_columns])
    assert list(no_be["Sequence"]) == ["AAAAA", "VVVVV"]
    assert list(no_be["Charge"]) == [2, 3]

    raw_diffs = dict(zip(no_be["Sequence"], no_be["pctD_diff"]))
    corrected_diffs = dict(zip(be["Sequence"], be["pctD_diff"]))
    assert raw_diffs == pytest.approx({"AAAAA": 20.0, "VVVVV": 10.0})
    assert corrected_diffs == pytest.approx({"AAAAA": 25.0, "VVVVV": 12.5})
    assert set(no_be["analysis_mode"]) == {"no_be"}
    assert set(no_be["normalization"]) == {"none"}
    assert set(be["analysis_mode"]) == {"be"}
    assert set(be["normalization"]) == {"standard"}

    # The workflow should request one render per comparison branch.
    assert len(pymol_calls) == 2
    for mode, command in zip(("no_be", "be"), pymol_calls):
        assert command[1:4] == ["-cq", "-r", str(Path(workflow.__file__).with_name("pymol_script.py"))]
        assert "--protein-only" in command
        assert command[command.index("--chain-filter") + 1] == "A"
        consensus_arg = Path(command[command.index("--value-col") - 1])
        assert consensus_arg == output / f"comparisons/A_vs_B/{mode}/tables/consensus_residue.csv"
        assert Path(command[command.index("--out-png") + 1]).parent == output / f"comparisons/A_vs_B/{mode}/structure"
