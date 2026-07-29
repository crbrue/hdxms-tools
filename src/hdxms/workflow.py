from __future__ import annotations

import math
import shutil
import hashlib
import platform
import sys
from importlib.metadata import PackageNotFoundError, version as package_version
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from .back_exchange import calculate_back_exchange_stats
from .consensus import consensus_from_csv
from .diffplots import plot_diff_bars, plot_diff_map, plot_aligned_diff
from .strip import export_consensus_full_manual_fixed_size, export_consensus_zoom_manual_fixed_size
from .quality import build_back_exchange_table, build_peptide_quality_table, quality_summary, write_quality_outputs
from .report import write_analysis_report
from .mass import peptide_monoisotopic_mass

PCT_CANDIDATES = ["Deut %", "% Max D", "%D", "Average %D", "Avg %D", "Percent D", "PctD", "Pct D"]
KEYS = ["Start", "End", "Sequence", "Charge"]


def _dependency_version(distribution: str, module_name: str | None = None) -> str:
    try:
        return package_version(distribution)
    except PackageNotFoundError:
        if module_name:
            module = __import__(module_name)
            return str(getattr(module, "__version__", "unknown"))
        return "unknown"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _pair_set(df: pd.DataFrame) -> set[tuple[str, Any]]:
    return set(zip(df["Sequence"].astype(str), df["Charge"]))


def _sequence_set(df: pd.DataFrame) -> set[str]:
    return set(df["Sequence"].astype(str))


def _valid_normalized_pairs(df: pd.DataFrame) -> set[tuple[str, Any]]:
    pc = _pct_col(df)
    valid = pd.to_numeric(df[pc], errors="coerce").replace([np.inf, -np.inf], np.nan).notna()
    valid &= pd.to_numeric(df["# Deut"], errors="coerce").replace([np.inf, -np.inf], np.nan).notna()
    return _pair_set(df.loc[valid])

def _pct_col(df: pd.DataFrame) -> str:
    for c in PCT_CANDIDATES:
        if c in df.columns:
            return c
    for c in df.columns:
        lc = str(c).lower()
        if "%d" in lc or ("deut" in lc and "%" in str(c)) or "percent" in lc:
            return c
    raise KeyError(f"Could not detect a percent-deuteration column. Columns: {list(df.columns)}")


def _clean(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "Confidence" in out.columns:
        out = out[out["Confidence"].astype(str).str.lower() != "low"]
    missing = [c for c in KEYS + ["# Deut"] if c not in out.columns]
    if missing:
        raise KeyError(f"HX Examiner export is missing required columns: {missing}")
    return out


def _read_clean(path: Path) -> pd.DataFrame:
    return _clean(pd.read_csv(path))


def _common_pairs(dfs: list[pd.DataFrame]) -> set[tuple[str, Any]]:
    return set.intersection(*[set(zip(df["Sequence"], df["Charge"])) for df in dfs])


def _filter_common(df: pd.DataFrame, common: set[tuple[str, Any]]) -> pd.DataFrame:
    mask = [(s, z) in common for s, z in zip(df["Sequence"], df["Charge"])]
    return df.loc[mask].sort_values(KEYS).reset_index(drop=True)


def _indexed(df: pd.DataFrame) -> pd.DataFrame:
    return df.set_index(["Sequence", "Charge"], drop=False)


def _normalize(df: pd.DataFrame, fd: pd.DataFrame, scale: float = 1.0) -> pd.DataFrame:
    """Normalize to FD. scale=fd_fraction/sample_fraction for empirical correction."""
    out = _indexed(df.copy())
    fd_i = _indexed(fd.copy())
    pc = _pct_col(out)
    fd_pc = _pct_col(fd_i)
    denominator = pd.to_numeric(fd_i[fd_pc], errors="coerce") * float(scale)
    out[pc] = pd.to_numeric(out[pc], errors="coerce") / denominator * 100.0
    # Preserve the original notebook behavior for the corrected # Deut field.
    out["# Deut"] = pd.to_numeric(out["# Deut"], errors="coerce") / denominator * 100.0
    return out.reset_index(drop=True).sort_values(KEYS)


def _condition_summary(frames: list[pd.DataFrame], output: Path) -> pd.DataFrame:
    # Inputs have already been restricted to the global common peptide/charge set.
    source = frames[0]
    optional = ["Search RT"] if "Search RT" in source.columns else []
    base = source[KEYS + optional].copy()
    peptide_ids = source["Peptide ID"] if "Peptide ID" in source.columns else pd.Series([None] * len(source), index=source.index)
    base["Peptide monoisotopic mass (uncharged)"] = [
        peptide_monoisotopic_mass(seq, charge, peptide_id)
        for seq, charge, peptide_id in zip(base["Sequence"], base["Charge"], peptide_ids)
    ]
    base["maxD"] = base["Sequence"].map(lambda s: max(len(str(s)) - 1 - str(s).count("P"), 0))
    reps = []
    for i, df in enumerate(frames, 1):
        col = f"Replicate {i}"
        base[col] = pd.to_numeric(df["# Deut"], errors="coerce").to_numpy()
        reps.append(col)
    base["Average #D"] = base[reps].mean(axis=1)
    base["SD"] = base[reps].std(axis=1)
    base.to_csv(output, index=False)
    return base


def _pct_stats(frames: list[pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for i, df in enumerate(frames, 1):
        pc = _pct_col(df)
        part = df[KEYS + [pc]].copy().rename(columns={pc: "pctD"})
        part["replicate"] = i
        rows.append(part)
    all_df = pd.concat(rows, ignore_index=True)
    return (all_df.groupby(KEYS, as_index=False)
            .agg(avg_pctD=("pctD", "mean"), sd_pctD=("pctD", "std"), n=("pctD", "count"))
            .sort_values(["Start", "End", "Charge"]).reset_index(drop=True))


def _resolve(root: Path, value: str) -> Path:
    p = Path(value)
    return p if p.is_absolute() else (root / p).resolve()


def run_manifest(manifest_path: str | Path) -> Path:
    manifest_path = Path(manifest_path).resolve()
    cfg = yaml.safe_load(manifest_path.read_text()) or {}
    root = _resolve(manifest_path.parent, str(cfg.get("base_dir", ".")))
    out = _resolve(root, str(cfg.get("output_dir", "hdxms_output")))
    out.mkdir(parents=True, exist_ok=True)
    provenance_dir = out / "provenance"
    provenance_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(manifest_path, provenance_dir / "manifest_used.yaml")

    sample_frac = float(cfg.get("sample_d2o_fraction", cfg.get("d2o_fraction", 1.0)))
    fd_frac = float(cfg.get("fd_d2o_fraction", 1.0))
    empirical = bool(cfg.get("apply_empirical_correction", False))
    if not (0 < sample_frac <= 1 and 0 < fd_frac <= 1):
        raise ValueError("sample_d2o_fraction and fd_d2o_fraction must be in (0, 1].")
    if empirical and math.isclose(sample_frac, fd_frac, rel_tol=0, abs_tol=1e-12):
        raise ValueError("apply_empirical_correction is true, but sample and FD D2O fractions are equal.")

    fd_path = _resolve(root, cfg["fully_deuterated"])
    datasets = cfg.get("datasets", cfg.get("conditions", {}))
    if len(datasets) < 1:
        raise ValueError("Manifest must define at least one dataset under 'datasets'.")
    dataset_titles = {name: str(spec.get("title", name)) for name, spec in datasets.items()}
    condition_paths = {name: [_resolve(root, f) for f in spec["files"]] for name, spec in datasets.items()}
    all_paths = [fd_path] + [p for paths in condition_paths.values() for p in paths]
    for p in all_paths:
        if not p.exists():
            raise FileNotFoundError(p)

    pipe = out / "pipeline_outputs"
    dirs = {name: pipe / name for name in [
        "step1_no_low", "step2_no_low_common_pair", "step3_fd_normalized",
        "step3_empirical_corrected", "step4_combined_corrected"
    ]}
    for d in dirs.values(): d.mkdir(parents=True, exist_ok=True)

    cleaned = {}
    for p in all_paths:
        df = _read_clean(p)
        cleaned[p] = df
        df.to_csv(dirs["step1_no_low"] / f"{p.stem}_no_low.csv", index=False)
    common = _common_pairs(list(cleaned.values()))

    # Exact raw-input inventory and immutable file fingerprints.
    inventory_rows = []
    path_roles: dict[Path, tuple[str, str, int | None]] = {fd_path: ("fully_deuterated", "FD", None)}
    for cond, paths in condition_paths.items():
        for rep_i, p in enumerate(paths, 1):
            path_roles[p] = ("sample", cond, rep_i)
    for p in all_paths:
        role, dataset_name, replicate = path_roles[p]
        st = p.stat()
        inventory_rows.append({
            "role": role,
            "dataset": dataset_name,
            "replicate": replicate,
            "input_name": p.name,
            "input_path": str(p),
            "sha256": _sha256(p),
            "size_bytes": st.st_size,
            "modified_utc": datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat(),
            "fd_control_name": fd_path.name if role == "sample" else "",
            "fd_control_path": str(fd_path) if role == "sample" else "",
            "fd_control_sha256": _sha256(fd_path) if role == "sample" else "",
            "sample_d2o_fraction": sample_frac if role == "sample" else "",
            "fd_d2o_fraction": fd_frac if role == "sample" else "",
        })
    pd.DataFrame(inventory_rows).to_csv(provenance_dir / "raw_input_inventory.csv", index=False)
    (provenance_dir / "software_environment.yaml").write_text(yaml.safe_dump({
        "python": sys.version,
        "platform": platform.platform(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "pyyaml": yaml.__version__,
        "pyteomics": _dependency_version("pyteomics", "pyteomics"),
    }, sort_keys=False))
    filtered = {}
    for p, df in cleaned.items():
        f = _filter_common(df, common)
        filtered[p] = f
        f.to_csv(dirs["step2_no_low_common_pair"] / f"{p.stem}_no_low_common_pair.csv", index=False)

    fd = filtered[fd_path]
    standard = {}
    empirical_frames = {}
    empirical_scale = fd_frac / sample_frac
    for p, df in filtered.items():
        std = _normalize(df, fd, scale=1.0)
        standard[p] = std
        std.to_csv(dirs["step3_fd_normalized"] / f"{p.stem}_fd_normalized.csv", index=False)
        if empirical:
            emp = _normalize(df, fd, scale=empirical_scale)
            empirical_frames[p] = emp
            emp.to_csv(dirs["step3_empirical_corrected"] / f"{p.stem}_empirical_corrected.csv", index=False)

    # Peptide attrition is reported independently for each input and dataset.
    attrition_rows = []
    global_common_sequences = set.intersection(*[_sequence_set(df) for df in cleaned.values()])
    for p in all_paths:
        role, dataset_name, replicate = path_roles[p]
        before_pairs = _pair_set(cleaned[p])
        after_common_pairs = _pair_set(filtered[p])
        after_be_pairs = _valid_normalized_pairs(standard[p])
        before_seq = _sequence_set(cleaned[p])
        after_common_seq = _sequence_set(filtered[p])
        after_be_seq = {seq for seq, _ in after_be_pairs}
        attrition_rows.append({
            "scope": "input_file",
            "role": role,
            "dataset": dataset_name,
            "replicate": replicate,
            "input_name": p.name,
            "peptide_charge_pairs_before_common_filter": len(before_pairs),
            "peptide_charge_pairs_after_common_filter": len(after_common_pairs),
            "peptide_charge_pairs_lost_common_filter": len(before_pairs - after_common_pairs),
            "peptide_charge_pairs_after_be_normalization": len(after_be_pairs),
            "peptide_charge_pairs_lost_during_be_normalization": len(after_common_pairs - after_be_pairs),
            "unique_sequences_before_common_filter": len(before_seq),
            "unique_sequences_after_common_filter": len(after_common_seq),
            "unique_sequences_lost_common_filter": len(before_seq - after_common_seq),
            "unique_sequences_after_be_normalization": len(after_be_seq),
            "unique_sequences_lost_during_be_normalization": len(after_common_seq - after_be_seq),
        })
    pd.DataFrame(attrition_rows).to_csv(provenance_dir / "peptide_attrition_by_input.csv", index=False)

    dataset_attrition = []
    for cond, paths in condition_paths.items():
        before_union = set.union(*[_pair_set(cleaned[p]) for p in paths])
        before_intersection = set.intersection(*[_pair_set(cleaned[p]) for p in paths])
        after_common = set.intersection(*[_pair_set(filtered[p]) for p in paths])
        after_be = set.intersection(*[_valid_normalized_pairs(standard[p]) for p in paths])
        dataset_attrition.append({
            "dataset": cond,
            "replicate_count": len(paths),
            "union_peptide_charge_pairs_before_filter": len(before_union),
            "within_dataset_common_pairs_before_global_filter": len(before_intersection),
            "global_common_pairs_retained": len(after_common),
            "pairs_lost_when_enforcing_global_common_set": len(before_intersection - after_common),
            "common_pairs_after_be_normalization": len(after_be),
            "pairs_lost_during_be_normalization": len(after_common - after_be),
            "fd_control_name": fd_path.name,
            "fd_control_sha256": _sha256(fd_path),
        })
    pd.DataFrame(dataset_attrition).to_csv(provenance_dir / "peptide_attrition_by_dataset.csv", index=False)

    selected = empirical_frames if empirical else standard
    selected_label = "empirical" if empirical else "standard"

    # Notebook-style combined replicate table, while keeping conditions explicit in column names.
    nonfd = [p for paths in condition_paths.values() for p in paths]
    base = filtered[nonfd[0]][KEYS].copy().set_index(["Sequence", "Charge"])
    for cond, paths in condition_paths.items():
        for i, p in enumerate(paths, 1):
            raw = _indexed(filtered[p]); std = _indexed(standard[p])
            base[f"# Deut_{cond}_rep{i}"] = raw["# Deut"]
            base[f"# Deut_std_{cond}_rep{i}"] = std["# Deut"]
            rpc = _pct_col(raw); spc = _pct_col(std)
            base[f"Deut %_{cond}_rep{i}"] = raw[rpc]
            base[f"Deut %_std_{cond}_rep{i}"] = std[spc]
            if empirical:
                emp = _indexed(empirical_frames[p]); epc = _pct_col(emp)
                base[f"# Deut_emp_corr_{cond}_rep{i}"] = emp["# Deut"]
                base[f"Deut %_emp_corr_{cond}_rep{i}"] = emp[epc]
    combined = base.reset_index().sort_values(KEYS)
    combined_path = dirs["step4_combined_corrected"] / "combined_replicates.csv"
    combined.to_csv(combined_path, index=False)

    summaries_dir = out / "condition_summaries"; summaries_dir.mkdir(exist_ok=True)
    stats = {}

    # v1.1 canonical dual exports. The BE and NO_BE tables use the exact same
    # _condition_summary schema; only the underlying numerical values differ.
    datasets_dir = out / "datasets"
    analysis_report_dir = out / "analysis_report"
    peptide_quality_dir = analysis_report_dir / "peptide_quality"
    back_exchange_dir = out / "back_exchange"
    datasets_dir.mkdir(exist_ok=True)
    analysis_report_dir.mkdir(exist_ok=True)
    back_exchange_dir.mkdir(exist_ok=True)

    fd_audit = build_back_exchange_table(fd, fd_frac)
    fd_audit.to_csv(back_exchange_dir / "peptide_back_exchange.csv", index=False)
    valid_be = pd.to_numeric(fd_audit["peptide_back_exchange_fraction"], errors="coerce").dropna()
    be_summary = {
        "n_peptides": int(len(fd_audit)),
        "n_valid": int(valid_be.size),
        "mean_back_exchange_fraction": float(valid_be.mean()) if len(valid_be) else None,
        "median_back_exchange_fraction": float(valid_be.median()) if len(valid_be) else None,
        "std_back_exchange_fraction": float(valid_be.std()) if len(valid_be) else None,
        "minimum_back_exchange_fraction": float(valid_be.min()) if len(valid_be) else None,
        "maximum_back_exchange_fraction": float(valid_be.max()) if len(valid_be) else None,
        "fd_control": fd_path.name,
        "fd_control_sha256": _sha256(fd_path),
        "fd_d2o_fraction": fd_frac,
    }
    (back_exchange_dir / "summary_statistics.yaml").write_text(yaml.safe_dump(be_summary, sort_keys=False))

    qc_tables = []
    dataset_report = {}
    qc_cfg = cfg.get("quality_control", {})
    for cond, paths in condition_paths.items():
        raw_frames = [filtered[p] for p in paths]
        be_frames = [standard[p] for p in paths]
        droot = datasets_dir / cond
        no_be_dir = droot / "no_be"
        be_dir = droot / "be"
        no_be_dir.mkdir(parents=True, exist_ok=True)
        be_dir.mkdir(parents=True, exist_ok=True)

        # MOST IMPORTANT FIRST SI TABLE contract: identical columns/order.
        no_be_summary = _condition_summary(raw_frames, no_be_dir / "SI_summary_table.csv")
        be_summary_table = _condition_summary(be_frames, be_dir / "SI_summary_table.csv")
        _condition_summary([selected[p] for p in paths], summaries_dir / f"{cond}_summary.csv")
        stats[cond] = _pct_stats([selected[p] for p in paths])

        if bool(qc_cfg.get("enabled", True)):
            qct = build_peptide_quality_table(cond, raw_frames, be_frames, fd_audit, qc_cfg)
            qc_tables.append(qct)
            dataset_report[cond] = {
                "title": dataset_titles[cond],
                "replicates": len(paths),
                "si_no_be": str((no_be_dir / "SI_summary_table.csv").relative_to(out)),
                "si_be": str((be_dir / "SI_summary_table.csv").relative_to(out)),
                "quality": quality_summary(qct),
            }
        else:
            dataset_report[cond] = {
                "title": dataset_titles[cond],
                "replicates": len(paths),
                "si_no_be": str((no_be_dir / "SI_summary_table.csv").relative_to(out)),
                "si_be": str((be_dir / "SI_summary_table.csv").relative_to(out)),
                "quality": {"enabled": False},
            }

    if qc_tables:
        all_qc = pd.concat(qc_tables, ignore_index=True)
        write_quality_outputs(all_qc, peptide_quality_dir)

    if empirical:
        be_dir = out / "back_exchange"; be_dir.mkdir(exist_ok=True)
        calculate_back_exchange_stats(combined_path, be_dir / "back_exchange_stats.csv", be_dir / "global_back_exchange_stats.csv")

    comparisons_dir = out / "comparisons"; comparisons_dir.mkdir(exist_ok=True)
    dpi = int(cfg.get("dpi", 600))
    for comp in cfg.get("comparisons", []):
        name = comp["name"]
        a = comp.get("dataset_a", comp.get("condition_a"))
        b = comp.get("dataset_b", comp.get("condition_b"))
        if a not in stats or b not in stats:
            raise KeyError(f"Comparison {name!r} references unknown datasets: {a!r}, {b!r}")
        cdir = comparisons_dir / name
        for sub in ["tables", "figures", "illustrator"]: (cdir / sub).mkdir(parents=True, exist_ok=True)
        sa = stats[a].rename(columns={"avg_pctD":f"avg_pctD_{a}","sd_pctD":f"sd_pctD_{a}","n":f"n_{a}"})
        sb = stats[b].rename(columns={"avg_pctD":f"avg_pctD_{b}","sd_pctD":f"sd_pctD_{b}","n":f"n_{b}"})
        merged = pd.merge(sa, sb, on=KEYS, how="inner").sort_values(["Start","End","Charge"])
        merged["pctD_diff"] = merged[f"avg_pctD_{a}"] - merged[f"avg_pctD_{b}"]
        pct_path = cdir / "tables" / "pctD_summary.csv"; merged.to_csv(pct_path, index=False)
        label = comp.get("label", f"{dataset_titles[a]} - {dataset_titles[b]}")
        plot_diff_bars(pct_path, cdir/"figures"/"pctD_diff_bars.png", label=label, dpi=dpi)
        plot_diff_map(pct_path, cdir/"figures"/"peptide_diff_map_stepwise.png", label=label, dpi=dpi)
        plot_diff_map(pct_path, cdir/"figures"/"peptide_diff_map_compact.png", label=label, compact=True, dpi=dpi)
        consensus_path = cdir/"tables"/"consensus_residue.csv"
        consensus_from_csv(pct_path, value_col="pctD_diff", offset=int(cfg.get("offset",0)), stat=cfg.get("stat","mean"), mincov=int(cfg.get("mincov",1)), smooth=int(cfg.get("smooth",0)), output=consensus_path)
        export_consensus_full_manual_fixed_size(consensus_path, cdir/"illustrator"/"consensus_full.svg")
        for z in cfg.get("zoom_regions", []):
            export_consensus_zoom_manual_fixed_size(consensus_path, int(z["start"]), int(z["end"]), cdir/"illustrator"/f"consensus_{z['name']}_{z['start']}-{z['end']}.svg")
        plot_aligned_diff(pct_path, cdir/"figures"/"consensus_domain_diffmap.png", offset=int(cfg.get("offset",0)), stat=cfg.get("stat","mean"), mincov=int(cfg.get("mincov",1)), smooth=int(cfg.get("smooth",0)), consensus_output=consensus_path, title=comp.get("title", name), dpi=dpi)

    (out / "run_metadata.yaml").write_text(yaml.safe_dump({
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": str(manifest_path),
        "manifest_sha256": _sha256(manifest_path),
        "output_directory": str(out),
        "normalization_used_for_primary_outputs": selected_label,
        "fully_deuterated_control": {
            "name": fd_path.name,
            "path": str(fd_path),
            "sha256": _sha256(fd_path),
        },
        "sample_d2o_fraction": sample_frac, "fd_d2o_fraction": fd_frac,
        "apply_empirical_correction": empirical,
        "global_common_peptide_charge_pairs": len(common),
        "global_common_unique_sequences": len(global_common_sequences),
        "datasets": {name: {
            "title": dataset_titles[name],
            "replicates": [p.name for p in condition_paths[name]],
            "replicate_paths": [str(p) for p in condition_paths[name]],
            "fd_control": fd_path.name,
        } for name in datasets},
        "provenance_files": [
            "provenance/raw_input_inventory.csv",
            "provenance/peptide_attrition_by_input.csv",
            "provenance/peptide_attrition_by_dataset.csv",
            "provenance/manifest_used.yaml",
            "provenance/software_environment.yaml",
        ],
    }, sort_keys=False))

    warnings = []
    for name, item in dataset_report.items():
        quality = item.get("quality", {})
        if quality.get("poor", 0):
            warnings.append({
                "level": "red",
                "dataset": name,
                "message": f"{quality['poor']} peptide/charge pairs scored Poor; review flagged_peptides.csv",
            })
        if quality.get("overcorrected_gt_110", 0):
            warnings.append({
                "level": "red",
                "dataset": name,
                "message": f"{quality['overcorrected_gt_110']} peptide/charge pairs exceeded 110% corrected occupancy",
            })

    inventory = [str(path.relative_to(out)) for path in out.rglob("*") if path.is_file()]
    write_analysis_report(
        analysis_report_dir / "analysis_report.yaml",
        manifest=str(manifest_path),
        normalization={
            "canonical_outputs": ["no_be", "be"],
            "sample_d2o_fraction": sample_frac,
            "fd_d2o_fraction": fd_frac,
            "primary_comparison_mode": selected_label,
            "empirical_d2o_correction_enabled": empirical,
            "si_table_schema_preserved_between_modes": True,
        },
        datasets=dataset_report,
        comparisons=cfg.get("comparisons", []),
        input_inventory_csv="provenance/raw_input_inventory.csv",
        peptide_attrition_csv="provenance/peptide_attrition_by_dataset.csv",
        output_inventory=inventory,
        warnings=warnings,
    )
    return out
