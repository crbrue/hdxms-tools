# HDX-MS Tools

# hdx-ms-tools

Installable HDX-MS utilities consolidated from the original HX Examiner, consensus-plotting, duplicate-removal, and PyMOL mapping notebooks/scripts.

## Install from GitHub

```bash
python -m pip install "hdx-ms-tools @ git+https://github.com/YOUR-USER/hdx-ms-tools.git"
```

For an editable development install:

```bash
git clone https://github.com/YOUR-USER/hdx-ms-tools.git
cd hdx-ms-tools
python -m pip install -e ".[dev,stats]"
```

PyMOL is invoked as an external executable and should be installed separately. The default executable name is `pymol`; override it with `hdxms color --pymol /path/to/pymol ...`.

## Consensus-based structure coloring

```bash
hdxms consensus Data_100s/pctD_summary/pctD_summary.csv \
  --value-col pctD_diff \
  --offset -23 \
  --mincov 2 \
  --smooth 5 \
  --resi-min 1 --resi-max 1014 \
  -o consensus_residue.csv

hdxms color PARP1_talazoparib_mono.pdb consensus_residue.csv \
  --value-col pctD_diff \
  --chain-filter A,B,C \
  --mincov 1 \
  --cartoon-trans 0.1 \
  --out-png PARP1_consensus_dHDX.png \
  --out-bfactor-pdb PARP1_consensus_dHDX_bfactor.pdb \
  --session PARP1_consensus_dHDX.pse
```

The consensus command expands each peptide over its residue span, combines overlapping peptide values, supports mean/median/length-weighted mean, applies minimum coverage, optionally smooths with an odd residue window, and explicitly retains uncovered residues as `NaN`.

## Illustrator-safe consensus strips

The consensus SVG commands use the manual exporters from
`consensus_plots_forIllustrator_02042026.ipynb`, including exact physical
sizes, one rectangle per residue, explicit white gaps, and
`preserveAspectRatio="none"`.

Full-length strip using the notebook defaults (7.4558 × 0.0751 inches):

```bash
hdxms strip-full consensus_residue.csv -o consensus_full.svg \
  --value-col pctD_diff
```

Zoom strip using the notebook defaults (1.6089 × 0.0392 inches):

```bash
hdxms strip-zoom consensus_residue.csv -o consensus_HD_zoom.svg \
  --value-col pctD_diff --resi-min 678 --resi-max 787
```

`hdxms strip` remains as a compatibility command. With no residue range it
calls the full Illustrator exporter; with both `--resi-min` and `--resi-max`
it calls the zoom Illustrator exporter.

## Remove duplicate HX Examiner peptide rows

```bash
hdxms deduplicate Data_100s/*.csv
```

For each `(Start, End, Sequence, Charge)` group, the row with the lowest non-missing `Search RT` is retained. Input files are changed atomically and receive a `.bak` backup by default.

## Back-exchange summaries

```bash
hdxms back-exchange combined_replicates.csv \
  --per-output back_exchange_per_peptide.csv \
  --global-output back_exchange_global.csv
```

## Color bins

- `< -30`: dark blue
- `-30 to < -15`: medium blue
- `-15 to < -5`: light blue
- `-5 to +5`: gray
- `> +5 to +15`: light red
- `> +15 to +30`: red
- `> +30`: dark red

## Repository layout

```text
src/hdxms/
  io.py             normalized CSV loading and residue offsets
  consensus.py      peptide-to-residue consensus calculation
  colors.py         one shared HDX color scale
  strip.py          Illustrator-safe SVG strips
  deduplicate.py    HX Examiner duplicate removal
  back_exchange.py  back-exchange summaries
  pymol_script.py   PyMOL renderer
  cli.py            hdxms command-line interface
```

## Differential HDX plots

Peptide/charge difference bars:

```bash
hdxms diff-bars pctD_summary.csv -o pctD_diff_bars.png \
  --value-col pctD_diff --label "Drug - Apo"
```

Stepwise peptide difference map (one row per peptide/charge):

```bash
hdxms diff-map pctD_summary.csv -o peptide_diff_map.png \
  --value-col pctD_diff --label "Drug - Apo" \
  --resi-min 1 --resi-max 1014
```

Use `--compact` to pack non-overlapping peptides into shared lanes. Zoom maps use the same command with a narrower `--resi-min/--resi-max` range.

Aligned consensus strip, PARP1 domain ruler, and compact peptide difference map:

```bash
hdxms diff-figure pctD_summary.csv -o aligned_diff_figure.png \
  --value-col pctD_diff --mincov 2 --smooth 5 \
  --resi-min 1 --resi-max 1014 \
  --consensus-output consensus_residue.csv \
  --title "PARP1 Δ%HDX — Drug vs Apo"
```

All differential visualizations use the same discrete blue/gray/red bins as the PyMOL structure coloring and Illustrator strip exports.

## End-to-end manifest workflow

The YAML manifest is the input to the complete workflow:

```bash
hdxms run examples/workflow_manifest.yaml
```

Standard fully-deuterated-control normalization is always generated. The optional
empirical D2O adjustment is controlled by:

```yaml
sample_d2o_fraction: 0.75
fd_d2o_fraction: 1.00
apply_empirical_correction: false
```

Set `apply_empirical_correction: true` only when the sample and FD-control D2O
fractions differ. The empirical denominator is scaled by
`fd_d2o_fraction / sample_d2o_fraction`. The workflow rejects an enabled
empirical correction when the two fractions are equal.

## Workflow specification documents

The repository includes the following normative documents:

- `docs/USAGE_AND_PROVENANCE.md` — installation, workflow use, and provenance requirements
- `docs/OUTPUT_CONTRACT.md` — canonical BE/NO_BE SI tables and required output tree
- `docs/MANIFEST_REFERENCE.md` — manifest fields and rules
- `docs/ANALYSIS_REPORT_SPEC.md` — run-level QC and provenance report
- `docs/PEPTIDE_QC_SPEC.md` — peptide scoring, flags, and red-flag export
- `docs/TESTING.md` — acceptance tests for validating the implementation

The original SI table remains the primary output. The BE and NO_BE versions must remain identical in structure, with all QC and provenance exported separately.

## Raw CSV filename helper

Generate an authoritative Python list of every `.csv` file in a raw-data directory:

```bash
python scripts/list_raw_csvs.py path/to/raw_data
```

This writes `raw_csv_files.py`:

```python
raw_files = [
    'APO_rep1.csv',
    'APO_rep2.csv',
    'APO_FD.csv',
]
```

To also request a **best-effort** manifest skeleton:

```bash
python scripts/list_raw_csvs.py path/to/raw_data \
    --suggest-manifest manifest_file_suggestion.yaml
```

The Python list is authoritative and contains every discovered CSV. Manifest grouping is deliberately conservative: ambiguous files and multiple FD candidates are left for manual assignment rather than silently guessed.

## Peptide monoisotopic masses and modifications

SI-table peptide masses are calculated with Pyteomics. The package preserves the
legacy workflow convention exactly: neutral monoisotopic peptide mass plus one
proton mass per observed charge state. Consequently, the historical column
`Peptide monoisotopic mass (uncharged)` is charge-state-specific and is not m/z.

When an HX Examiner export includes a `Peptide ID` containing valid ProForma
modification notation, that modified peptidoform is used for the mass calculation.
The unmodified sequence obtained from the ProForma identifier must match the
`Sequence` column; otherwise the workflow stops with an explicit error. Plain or
non-ProForma peptide identifiers do not override the `Sequence` column.
