#!/usr/bin/env bash
# =============================================================================
# Project : kelpPR (Dongshan Bay kelp aquaculture metagenome)
# Author  : Fei Zhang
# Date    : 2026-09-14
# Version : 1.0
# Input   : 02_mag/bin_refinement/metawrap_70_10_bins/*.fa (131 MAGs)
#           00_data/geneset/PF01036.hmm (Pfam 38.2 proteorhodopsin HMM)
# Output  : 05_mag_pr/mag_proteins/*.faa        (prodigal per MAG)
#           05_mag_pr/mag_pr_hits.tbl           (hmmsearch --tblout, E<=1e-10)
#           05_mag_pr/mag_pr_summary.tsv        (MAG, n_PR, PR_ids, taxonomy)
# Notes   : independent evidence for PR hosts, complements step7 read-mapping.
#           Length filter 150-450 aa mirrors step2 of the PR pipeline.
# Resources: prodigal single-thread x131 (fast); hmmsearch --cpu 8, nice -n 10
# =============================================================================
set -euo pipefail

PROJ="${PROJECT_DIR:-$HOME/projects/kelpPR}"
OUT=$PROJ/05_mag_pr
BINS=$PROJ/02_mag/bin_refinement/metawrap_70_10_bins
HMM=$PROJ/00_data/geneset/PF01036.hmm
mkdir -p "$OUT/mag_proteins"
cd "$OUT"

echo "[$(date '+%F %T')] prodigal on 131 MAGs"
for b in "$BINS"/*.fa; do
  bin=$(basename "$b" .fa)
  [ -s "mag_proteins/${bin}.faa" ] && continue
  nice -n 10 prodigal -i "$b" -a "mag_proteins/${bin}.faa" -p meta -q
done
echo "[$(date '+%F %T')] prodigal done: $(ls mag_proteins/*.faa | wc -l) proteomes"

cat mag_proteins/*.faa > all_mag_proteins.faa
echo "[$(date '+%F %T')] total proteins: $(grep -c '^>' all_mag_proteins.faa)"

nice -n 10 hmmsearch --cpu 8 -E 1e-10 --tblout mag_pr_hits.tbl "$HMM" all_mag_proteins.faa > /dev/null
echo "[$(date '+%F %T')] hmmsearch done"

python3 step9_summarize.py
echo "[$(date '+%F %T')] Step9 done -> mag_pr_summary.tsv"
