from pathlib import Path

import pytest

import hdxms.cli as cli


def test_color_cli_forwards_all_structure_rendering_options(monkeypatch, tmp_path: Path):
    structure = tmp_path / "protein.pdb"
    values = tmp_path / "consensus.csv"
    png = tmp_path / "colored.png"
    session = tmp_path / "colored.pse"
    bfactor = tmp_path / "colored_bfactor.pdb"
    structure.write_text("END\n")
    values.write_text("resi,pctD_diff\n1,2.5\n")

    captured = {}

    def fake_call(command):
        captured["command"] = [str(value) for value in command]
        return 0

    monkeypatch.setattr(cli.subprocess, "call", fake_call)

    with pytest.raises(SystemExit) as exc:
        cli.main([
            "color",
            str(structure),
            str(values),
            "--pymol", "pymol-custom",
            "--value-col", "custom_value",
            "--chain-filter", "A,B",
            "--offset", "5",
            "--offset-per-chain", "A:10,B:-2",
            "--mincov", "3",
            "--cartoon-trans", "0.25",
            "--protein-only",
            "--out-png", str(png),
            "--session", str(session),
            "--out-bfactor-pdb", str(bfactor),
        ])

    assert exc.value.code == 0
    command = captured["command"]
    assert command[:3] == ["pymol-custom", "-cq", "-r"]
    assert command[3] == str(Path(cli.__file__).with_name("pymol_script.py"))
    assert command[4:7] == ["--", str(structure), str(values)]

    expected_pairs = {
        "--value-col": "custom_value",
        "--chain-filter": "A,B",
        "--offset": "5",
        "--offset-per-chain": "A:10,B:-2",
        "--mincov": "3",
        "--cartoon-trans": "0.25",
        "--out-png": str(png),
        "--session": str(session),
        "--out-bfactor-pdb": str(bfactor),
    }
    for flag, expected in expected_pairs.items():
        assert command[command.index(flag) + 1] == expected
    assert "--protein-only" in command


def test_color_cli_returns_pymol_exit_code(monkeypatch, tmp_path: Path):
    structure = tmp_path / "protein.pdb"
    values = tmp_path / "consensus.csv"
    structure.write_text("END\n")
    values.write_text("resi,pctD_diff\n1,0.0\n")
    monkeypatch.setattr(cli.subprocess, "call", lambda command: 7)

    with pytest.raises(SystemExit) as exc:
        cli.main(["color", str(structure), str(values)])

    assert exc.value.code == 7
