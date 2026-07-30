from __future__ import annotations

import importlib
import sys
import types
from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest


class _Atom:
    def __init__(self, chain: str, resi: str):
        self.chain = chain
        self.resi = resi


class _Model:
    def __init__(self, atoms):
        self.atom = atoms


class FakeCmd:
    def __init__(self):
        self.calls = []
        self.saved = []
        self.model = _Model([
            _Atom("A", "1"),
            _Atom("A", "2"),
            _Atom("A", "3"),
            _Atom("B", "11"),
            _Atom("B", "12"),
        ])

    def __getattr__(self, name):
        def method(*args, **kwargs):
            self.calls.append((name, args, kwargs))
            if name == "get_model":
                return self.model
            if name == "save":
                self.saved.append(args)
        return method


class FakeUtil:
    def __init__(self):
        self.calls = []

    def cbaw(self, selection):
        self.calls.append(selection)


def _import_pymol_script(monkeypatch):
    fake_cmd = FakeCmd()
    fake_util = FakeUtil()
    pymol_module = types.ModuleType("pymol")
    pymol_module.cmd = fake_cmd
    pymol_module.util = fake_util
    monkeypatch.setitem(sys.modules, "pymol", pymol_module)
    sys.modules.pop("hdxms.pymol_script", None)
    module = importlib.import_module("hdxms.pymol_script")
    return module, fake_cmd, fake_util


def _args(tmp_path: Path, csv_path: Path, **overrides):
    defaults = {
        "structure": str(tmp_path / "protein.pdb"),
        "csv": str(csv_path),
        "value_col": "pctD_diff",
        "chain_filter": None,
        "offset": 0,
        "offset_per_chain": None,
        "mincov": 1,
        "cartoon_trans": 0.0,
        "protein_only": True,
        "out_png": str(tmp_path / "mapped.png"),
        "session": str(tmp_path / "mapped.pse"),
        "out_bfactor_pdb": str(tmp_path / "mapped_bfactor.pdb"),
        "width": 800,
        "height": 600,
        "dpi": 150,
    }
    defaults.update(overrides)
    Path(defaults["structure"]).write_text("END\n")
    return Namespace(**defaults)


def test_pymol_mapping_writes_expected_chain_specific_bfactors(monkeypatch, tmp_path: Path):
    module, fake_cmd, _ = _import_pymol_script(monkeypatch)
    csv_path = tmp_path / "values.csv"
    pd.DataFrame([
        {"Start": 1, "End": 2, "Sequence": "AA", "Charge": 2, "chain": "A", "pctD_diff": 4.0},
        {"Start": 11, "End": 12, "Sequence": "BB", "Charge": 2, "chain": "B", "pctD_diff": -3.0},
    ]).to_csv(csv_path, index=False)

    module.run(_args(tmp_path, csv_path, chain_filter="A,B"))

    alter_calls = [args for name, args, _ in fake_cmd.calls if name == "alter"]
    # The script intentionally resets all protein B-factors before assigning
    # mapped residue values, preventing stale structure values from surviving.
    assert ("protein and polymer.protein", "b=0.0") in alter_calls
    assert any(args[0].endswith("and chain A and resi 1") and args[1] == "b=4.0" for args in alter_calls)
    assert any(args[0].endswith("and chain A and resi 2") and args[1] == "b=4.0" for args in alter_calls)
    assert any(args[0].endswith("and chain B and resi 11") and args[1] == "b=-3.0" for args in alter_calls)
    assert any(args[0].endswith("and chain B and resi 12") and args[1] == "b=-3.0" for args in alter_calls)
    assert not any("resi 3" in args[0] and args[1] != "b=0.0" for args in alter_calls)
    assert any(name == "save" and args[0].endswith("mapped_bfactor.pdb") for name, args, _ in fake_cmd.calls)


def test_pymol_mapping_applies_chain_offsets_and_minimum_coverage(monkeypatch, tmp_path: Path):
    module, fake_cmd, _ = _import_pymol_script(monkeypatch)
    csv_path = tmp_path / "values.csv"
    pd.DataFrame([
        {"Start": 1, "End": 2, "Sequence": "AA", "Charge": 2, "chain": "A", "pctD_diff": 2.0},
        {"Start": 1, "End": 1, "Sequence": "A", "Charge": 3, "chain": "A", "pctD_diff": 6.0},
        {"Start": 1, "End": 2, "Sequence": "BB", "Charge": 2, "chain": "B", "pctD_diff": -4.0},
    ]).to_csv(csv_path, index=False)

    module.run(_args(
        tmp_path,
        csv_path,
        chain_filter="A,B",
        offset_per_chain="A:0,B:10",
        mincov=2,
    ))

    alter_calls = [args for name, args, _ in fake_cmd.calls if name == "alter"]
    # A1 has two-peptide coverage and receives their mean; A2 and B residues do not meet mincov=2.
    assert any(args[0].endswith("and chain A and resi 1") and args[1] == "b=4.0" for args in alter_calls)
    assert not any("chain A and resi 2" in args[0] and args[1] != "b=0.0" for args in alter_calls)
    assert not any("chain B and resi" in args[0] and args[1] != "b=0.0" for args in alter_calls)


def test_ranges_collapses_contiguous_residues(monkeypatch):
    module, _, _ = _import_pymol_script(monkeypatch)
    assert module._ranges([5, 1, 2, 3, 7, 7]) == ["1-3", "5", "7"]
    assert module._ranges([]) == []
