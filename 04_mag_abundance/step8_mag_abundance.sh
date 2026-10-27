#!/usr/bin/env bash
# =============================================================================
# Project : kelpPR (Dongshan Bay kelp aquaculture metagenome)
# Author  : Fei Zhang
# Date    : 2026-09-14
# Version : 1.0
# Input   : 02_mag/bin_refinement/metawrap_70_10_bins/*.fa (131 MAGs)
#           02_mag/coassembly_megahit/final.contigs.fa     (co-assembly)
#           00_data/fastq_metawrap/S01_R000{1..6}_{1,2}.fastq (softlinks to clean reads)
# Output  : 04_mag_abundance/quant/bin_abundance_table.tab (salmon TPM per MAG x sample)
#           04_mag_abundance/mag_hd_vs_ck.tsv              (group means, HD/CK ratio)
# Notes   : metaWRAP quant_bins (salmon). R0001-3 = CK, R0004-6 = HD (verify!).
#           n=3 vs 3: report descriptive statistics + Mann-Whitney as limit.
# Resources: salmon -t 8 per run (sequential samples), nice -n 10
# =============================================================================
set -euo pipefail

PROJ="${PROJECT_DIR:-$HOME/projects/kelpPR}"
OUT=$PROJ/04_mag_abundance
BINS=$PROJ/02_mag/bin_refinement/metawrap_70_10_bins
ASM=$PROJ/02_mag/coassembly_megahit/final.contigs.fa
FQ=$PROJ/00_data/fastq_metawrap
mkdir -p "$OUT"
cd "$OUT"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate metawrap-mg

READS=""
for i in 1 2 3 4 5 6; do
  READS="$READS $FQ/S01_R000${i}_1.fastq $FQ/S01_R000${i}_2.fastq"
done

if [ ! -s quant/bin_abundance_table.tab ]; then
  echo "[$(date '+%F %T')] metaWRAP quant_bins start"
  nice -n 10 metawrap quant_bins -b "$BINS" -o quant -a "$ASM" -t 8 $READS
fi
conda deactivate
echo "[$(date '+%F %T')] quant_bins done"

python3 step8_compare.py
echo "[$(date '+%F %T')] Step8 done -> mag_hd_vs_ck.tsv"
