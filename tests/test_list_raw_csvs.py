from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

SCRIPT = Path(__file__).parents[1] / "scripts" / "list_raw_csvs.py"
spec = spec_from_file_location("list_raw_csvs", SCRIPT)
module = module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(module)


def test_csv_names_and_python_output(tmp_path):
    (tmp_path / "B.csv").write_text("x\n")
    (tmp_path / "a.csv").write_text("x\n")
    (tmp_path / "ignore.txt").write_text("x\n")
    names = module.csv_names(tmp_path)
    assert names == ["a.csv", "B.csv"]
    assert module.python_list_text(names) == "raw_files = [\n    'a.csv',\n    'B.csv',\n]\n"


def test_manifest_inference_is_conservative():
    names = ["APO_rep1.csv", "APO_rep2.csv", "APO_FD.csv", "APO_fullD.csv"]
    text = module.yaml_suggestion(names)
    assert "Multiple FD candidates detected" in text
    assert "Files requiring manual assignment" in text
    assert "fd_file:" not in text
