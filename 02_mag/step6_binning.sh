#!/bin/bash
# =============================================================================
# step6_binning.sh - MAG binning relay pipeline
# Project : kelpPR metagenomes (dataset 1: Dongshan Bay HD vs CK, n = 3 vs 3)
# Author  : Fei Zhang
# Date    : 2026-09-11 | Version: 1.2 (fix: the metaWRAP source requires names ending strictly in _1.fastq/_2.fastq; a .gz suffix is rejected)
# Input   : 02_mag/coassembly_megahit/final.contigs.fa  (MEGAHIT co-assembly, --min-contig-len 1000)
#           00_data/fastq/Unknown_CM794-S01-001R000{1..6}_good_{1,2}.fq.gz  (6 样本 clean reads)
# Output  : 02_mag/binning/            (metaBAT2 + MaxBin2 + CONCOCT 原始 bins)
#           02_mag/bin_refinement/     (three-binner refinement, completeness >= 70, contamination < 10)
#           02_mag/checkm/qa.tsv       (CheckM lineage_wf quality table)
#           02_mag/gtdbtk/             (GTDB-Tk v2.5.2 + r226 classification)
# Resources: binning/refinement -t 16; checkm -t 8; gtdbtk --cpus 8 (+ --scratch_dir to limit memory);
#            nice -n 10 throughout (shared-server policy, total <= 32 cores)
# Usage   : nohup bash step6_binning.sh > ../99_logs/step6_pipeline.log 2>&1 &
# Notes   : resumable - outputs that already exist are skipped; each stage validates its result and aborts if it is missing
# =============================================================================
set -uo pipefail   # no -e: metaWRAP sub-steps may fail partially, so each stage is validated explicitly

PROJ="${PROJECT_DIR:-$HOME/projects/kelpPR}"
ASM=$PROJ/02_mag/coassembly_megahit/final.contigs.fa
FQ=$PROJ/00_data/fastq_metawrap  # metaWRAP requires *_1.fastq / *_2.fastq naming; the original .fq.gz was rejected (fixed 2026-09-11)
OUT=$PROJ/02_mag
LOG=$PROJ/99_logs
T_BIN=16   # binning / refinement threads
T_QA=8     # checkm / gtdbtk threads
mkdir -p $LOG

ts() { date '+%F %T'; }

# ---------- Stage 0: prerequisites ----------
if [ ! -f "$ASM" ]; then
  echo "[$(ts)] [abort] co-assembly result not found: $ASM (MEGAHIT not finished)"
  exit 1
fi
echo "[$(ts)] co-assembly ready: $(grep -c '>' $ASM) contigs"

# read pairing (metaWRAP expects _1 and _2 adjacent for each sample)
READS=""
for i in 1 2 3 4 5 6; do
  r1=$FQ/S01_R000${i}_1.fastq  # symlinks point to .fq.gz (bwa reads gz transparently)
  r2=$FQ/S01_R000${i}_2.fastq
  [ -f "$r1" ] && [ -f "$r2" ] || { echo "[$(ts)] [abort] missing fastq: R000${i}"; exit 1; }
  READS="$READS $r1 $r2"
done

source ~/miniconda3/etc/profile.d/conda.sh

# ---------- Stage 1: three binners (metaBAT2 + MaxBin2 + CONCOCT) ----------
if [ -d "$OUT/binning/metabat2_bins" ] && [ "$(ls $OUT/binning/metabat2_bins/*.fa 2>/dev/null | wc -l)" -gt 0 ]; then
  echo "[$(ts)] [skip] Stage1 three-binner binning already done"
else
  echo "[$(ts)] Stage1: metaWRAP three-binner binning started (-t $T_BIN, expect several hours)"
  conda activate metawrap-mg
  nice -n 10 metawrap binning -o $OUT/binning -t $T_BIN -m 250 -a "$ASM" \
    --metabat2 --maxbin2 --concoct --universal $READS > $LOG/step6_binning.log 2>&1
  n_mb=$(ls $OUT/binning/metabat2_bins/*.fa 2>/dev/null | wc -l)
  n_mx=$(ls $OUT/binning/maxbin2_bins/*.fa 2>/dev/null | wc -l)
  n_cc=$(ls $OUT/binning/concoct_bins/*.fa 2>/dev/null | wc -l)
  echo "[$(ts)] Stage1 done: metabat2=$n_mb maxbin2=$n_mx concoct=$n_cc bins"
  if [ "$n_mb" -eq 0 ] && [ "$n_mx" -eq 0 ] && [ "$n_cc" -eq 0 ]; then
    echo "[$(ts)] [abort] all three binners produced nothing, see $LOG/step6_binning.log"
    exit 1
  fi
fi

# ---------- Stage 2: bin_refinement（-c 70 -x 10）----------
REFINED=$OUT/bin_refinement/metawrap_70_10_bins
if [ -d "$REFINED" ] && [ "$(ls $REFINED/*.fa 2>/dev/null | wc -l)" -gt 0 ]; then
  echo "[$(ts)] [skip] Stage2 bin_refinement already done"
else
  echo "[$(ts)] Stage2: bin_refinement started (completeness >= 70, contamination < 10)"
  conda activate metawrap-mg 2>/dev/null
  ARGS=""
  [ -d "$OUT/binning/metabat2_bins" ] && ARGS="$ARGS -A $OUT/binning/metabat2_bins"
  [ -d "$OUT/binning/maxbin2_bins" ] && ARGS="$ARGS -B $OUT/binning/maxbin2_bins"
  [ -d "$OUT/binning/concoct_bins" ] && ARGS="$ARGS -C $OUT/binning/concoct_bins"
  nice -n 10 metawrap bin_refinement -o $OUT/bin_refinement -t $T_BIN $ARGS \
    -c 70 -x 10 > $LOG/step6_refinement.log 2>&1
  n_ref=$(ls $REFINED/*.fa 2>/dev/null | wc -l)
  echo "[$(ts)] Stage2 done: $n_ref bins after refinement (70/10 threshold)"
  if [ "$n_ref" -eq 0 ]; then
    echo "[$(ts)] [warning] no bins at the 70/10 threshold - consider relaxing to 50/10 and re-running this stage, see $LOG/step6_refinement.log"
    exit 2
  fi
fi

# ---------- Stage 3: CheckM quality check (dedicated checkm environment) ----------
if [ -s "$OUT/checkm/qa.tsv" ]; then
  echo "[$(ts)] [skip] Stage3 CheckM already done"
else
  echo "[$(ts)] Stage3: CheckM lineage_wf started (-t $T_QA)"
  conda activate checkm
  mkdir -p $OUT/checkm
  nice -n 10 checkm lineage_wf -t $T_QA -x fa "$REFINED" $OUT/checkm \
    --tab_table -f $OUT/checkm/qa.tsv > $LOG/step6_checkm.log 2>&1
  if [ ! -s "$OUT/checkm/qa.tsv" ]; then
    echo "[$(ts)] [abort] CheckM produced no qa.tsv, see $LOG/step6_checkm.log"
    exit 1
  fi
  echo "[$(ts)] Stage3 done: $(($(wc -l < $OUT/checkm/qa.tsv) - 1)) bins quality-checked"
fi

# ---------- Stage 4: GTDB-Tk classification (v2.5.2 + r226, scratch mode to limit memory) ----------
if [ -s "$OUT/gtdbtk/gtdbtk.bac120.summary.tsv" ]; then
  echo "[$(ts)] [skip] Stage4 GTDB-Tk already done"
else
  echo "[$(ts)] Stage4: GTDB-Tk classify_wf started (--cpus $T_QA, --scratch_dir to limit memory)"
  conda deactivate 2>/dev/null
  source $PROJ/gtdbtk_env.sh
  mkdir -p $OUT/gtdbtk/scratch
  nice -n 10 gtdbtk classify_wf --genome_dir "$REFINED" --out_dir $OUT/gtdbtk \
    --cpus $T_QA -x fa --scratch_dir $OUT/gtdbtk/scratch > $LOG/step6_gtdbtk.log 2>&1
  if [ ! -s "$OUT/gtdbtk/gtdbtk.bac120.summary.tsv" ]; then
    echo "[$(ts)] [abort] GTDB-Tk produced no summary, see $LOG/step6_gtdbtk.log"
    exit 1
  fi
  echo "[$(ts)] Stage4 done: $(($(wc -l < $OUT/gtdbtk/gtdbtk.bac120.summary.tsv) - 1)) bacterial bins classified"
fi

echo "[$(ts)] ===== Step6 binning relay complete ====="
