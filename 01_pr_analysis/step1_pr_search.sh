#!/bin/bash
# =====================================================================
# PR identification pipeline - dataset 1 (kelp-farming vs control metagenomes)
# Step 1: hmmsearch PF01036 (Pfam 38.2) against the non-redundant gene catalogue
# Step 2: extract hits (e-value <= 1e-10, coverage >= 70%)
# Step 3: summary statistics
# Environment: conda activate Metagenomic (HMMER 3.4)
# =====================================================================
set -euo pipefail

BASE="${PROJECT_DIR:-$HOME/projects/kelpPR}"
DATA=$BASE/00_data/geneset
OUT=$BASE/01_pr_analysis
LOG=$BASE/99_logs

mkdir -p $OUT $LOG

echo "[$(date)] Step 1: hmmsearch PF01036 vs gene catalogue (2,482,502 proteins)..."
# Multi-threading happens inside a single hmmsearch call (--cpu capped for the shared server)
hmmsearch --cpu 8 \
  -E 1e-10 \
  --domT 0 \
  --incE 1e-10 \
  -o $LOG/hmmsearch_pf01036.out \
  --tblout $OUT/hmmsearch_pf01036.tbl \
  $DATA/PF01036.hmm \
  $DATA/geneset_protein.fasta

echo "[$(date)] Total hits (domtbl not needed, using tbl):"
grep -v "^#" $OUT/hmmsearch_pf01036.tbl | grep -c . || true

echo "[$(date)] Step 2: extract hit sequences + coverage filter..."
# Extract hit gene IDs from the tblout (column 1 = target name)
grep -v "^#" $OUT/hmmsearch_pf01036.tbl | awk '$5 <= 1e-10 {print $1}' | sort -u > $OUT/pr_hit_ids.txt
echo "Hits passing E<=1e-10: $(wc -l < $OUT/pr_hit_ids.txt)"

# Extract sequences with seqkit
seqkit grep -f $OUT/pr_hit_ids.txt $DATA/geneset_protein.fasta > $OUT/pr_candidates.fasta

# Length filter: keep 150-450 aa (typical 7-TM rhodopsin length range)
seqkit seq -m 150 -M 450 $OUT/pr_candidates.fasta > $OUT/pr_candidates_len.fasta

echo "[$(date)] Candidate PR sequences (150-450 aa): $(grep -c '^>' $OUT/pr_candidates_len.fasta)"

echo "[$(date)] Step 3: summary"
echo "=== PR search: preliminary results ===" > $OUT/pr_search_summary.txt
echo "Total sequences in the gene catalogue: 2482502" >> $OUT/pr_search_summary.txt
echo "PF01036 HMM hits (E<=1e-10): $(wc -l < $OUT/pr_hit_ids.txt)" >> $OUT/pr_search_summary.txt
echo "Candidates after the length filter (150-450 aa): $(grep -c '^>' $OUT/pr_candidates_len.fasta)" >> $OUT/pr_search_summary.txt

cat $OUT/pr_search_summary.txt
echo "[$(date)] Done. Next: TMHMM/Phobius 7-TM filtering + Lys231 check + CD-HIT OTU clustering"
