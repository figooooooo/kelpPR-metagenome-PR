#!/usr/bin/env bash
# =============================================================================
# Project : kelpPR (Dongshan Bay kelp aquaculture metagenome)
# Author  : Fei Zhang
# Date    : 2026-09-14
# Version : 1.0
# Input   : 00_data/geneset/geneset_gene.fasta      (non-redundant gene catalogue, nt)
#           01_pr_analysis/pr_filtered.fasta        (139 PR protein records, IDs)
#           02_mag/coassembly_megahit/final.contigs.fa (MEGAHIT co-assembly)
#           02_mag/bin_refinement/metawrap_70_10_bins/*.fa (131 MAGs)
#           02_mag/gtdbtk/classify/gtdbtk.*.summary.tsv  (GTDB-Tk r226 taxonomy)
# Output  : 03_pr_mag/pr_139_nt.fasta               (139 PR nucleotide sequences)
#           03_pr_mag/pr_vs_coassembly.blast6       (blastn megablast hits)
#           03_pr_mag/contig2bin.tsv                (co-assembly contig -> MAG)
#           03_pr_mag/pr_mag_mapping.tsv            (PR -> contig -> MAG -> taxonomy)
# Notes   : PR genes derive from Biomarker per-sample gene catalogue (ctg=CK-x|ctg*),
#           MAGs from MEGAHIT co-assembly (k141_*); sequence-level mapping required.
#           Thresholds: megablast, pident>=95, qcovs>=80, best hit by bitscore.
# Resources: blastn -num_threads 8, nice -n 10 (shared server policy)
# =============================================================================
set -euo pipefail

PROJ="${PROJECT_DIR:-$HOME/projects/kelpPR}"
OUT=$PROJ/03_pr_mag
ASM=$PROJ/02_mag/coassembly_megahit/final.contigs.fa
BINS=$PROJ/02_mag/bin_refinement/metawrap_70_10_bins
GTDB=$PROJ/02_mag/gtdbtk/classify
mkdir -p "$OUT"
cd "$OUT"

echo "[$(date '+%F %T')] Step7 PR->MAG mapping start"

# --- 1. Extract 139 PR nucleotide sequences ---------------------------------
if [ ! -s pr_139_nt.fasta ]; then
  grep '^>' "$PROJ/01_pr_analysis/pr_filtered.fasta" | sed 's/^>//; s/ .*//' | sort -u > pr_139_ids.txt
  seqkit grep -f pr_139_ids.txt "$PROJ/00_data/geneset/geneset_gene.fasta" > pr_139_nt.fasta
fi
echo "[$(date '+%F %T')] PR nt sequences: $(grep -c '^>' pr_139_nt.fasta)"

# --- 2. BLAST database of co-assembly ---------------------------------------
if [ ! -s coassembly_db.ndb ] && [ ! -s coassembly_db.nsq ]; then
  nice -n 10 makeblastdb -in "$ASM" -dbtype nucl -out coassembly_db -parse_seqids
fi
echo "[$(date '+%F %T')] BLAST db ready"

# --- 3. blastn PR nt vs co-assembly ------------------------------------------
nice -n 10 blastn -query pr_139_nt.fasta -db coassembly_db \
  -task megablast -perc_identity 95 \
  -num_threads 8 -max_target_seqs 5 \
  -outfmt '6 qseqid sseqid pident length qlen slen qcovs evalue bitscore' \
  -out pr_vs_coassembly.blast6
echo "[$(date '+%F %T')] blastn done: $(wc -l < pr_vs_coassembly.blast6) hits"

# --- 4. contig -> bin map -----------------------------------------------------
> contig2bin.tsv
for b in "$BINS"/*.fa; do
  bin=$(basename "$b" .fa)
  grep '^>' "$b" | sed 's/^>//; s/ .*//' | awk -v B="$bin" '{print $1"\t"B}'
done > contig2bin.tsv
echo "[$(date '+%F %T')] contig2bin: $(wc -l < contig2bin.tsv) contigs in MAGs"

# --- 5. Join PR -> contig -> MAG -> taxonomy ---------------------------------
python3 step7_join.py
echo "[$(date '+%F %T')] Step7 done -> pr_mag_mapping.tsv"
