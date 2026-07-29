#!/usr/bin/env bash
set -euo pipefail

hdxms consensus pctD_summary.csv \
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
