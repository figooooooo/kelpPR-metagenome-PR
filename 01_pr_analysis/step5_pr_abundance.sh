#!/bin/bash
# Step 5: PR abundance normalisation (GPM method)
# Prerequisites:
#   1) ~/projects/kelpPR/00_data/geneset/geneset_gene.fasta uploaded and decompressed
#   2) ~/projects/kelpPR/00_data/fastq/*.fq.gz transferred for all 6 samples
#   3) ~/projects/kelpPR/01_pr_analysis/pr_pass_ids.txt (139 PR gene IDs)
#   4) ~/projects/kelpPR/00_data/geneset/scg_ids.txt (recA/rpoB/rplB gene IDs)
#   5) ~/projects/kelpPR/00_data/geneset/ko_scg_map.tsv (gene ID -> KO mapping, used for classification)
# Usage : bash step5_pr_abundance.sh  (on the server: nohup bash step5_pr_abundance.sh &)
# Resource budget (agreed 2026-09-10):
#   bowtie2 keeps 8 threads per task; 6 samples run 3 at a time -> 24 cores peak (~19% of the 128-core server)
#   every compute process runs under nice -n 10 so other users are not starved; samtools view drops to 4 threads (compression is I/O bound)
set -euo pipefail

# Explicit PATH: Metagenomic environment (samtools/bowtie2) plus system tools
# Activate the required conda environment before running, e.g.
#   conda activate Metagenomic

PROJ="${PROJECT_DIR:-$HOME/projects/kelpPR}"
GENESET=$PROJ/00_data/geneset/geneset_gene.fasta
FASTQ_DIR=$PROJ/00_data/fastq
OUT=$PROJ/03_abundance
LOG=$PROJ/99_logs
THREADS=8           # threads per alignment task
SAM_THREADS=4       # samtools view threads (compression is I/O bound; 4 is enough)
CONC=3              # samples aligned concurrently (3 x 8 = 24 cores peak)
mkdir -p $OUT $OUT/bam $LOG

# Prerequisite 0: decompress the nucleotide gene catalogue if only .gz is present
if [ ! -f "$GENESET" ] && [ -f "${GENESET}.gz" ]; then
  echo "[$(date '+%F %T')] Decompressing the nucleotide gene catalogue"
  gunzip -k "${GENESET}.gz"
fi

echo "[$(date '+%F %T')] Step5: building the bowtie2 index"
# bowtie2 index (one-off, roughly 10-20 min)
if [ ! -f "$GENESET.bti" ] && [ ! -f "${GENESET}.1.bt2l" ] && [ ! -f "${GENESET}.1.bt2" ]; then
  nice -n 10 bowtie2-build --threads $THREADS "$GENESET" "$GENESET" > $LOG/step5_bt2build.log 2>&1
fi

# Align one sample (run under nice -n 10)
map_one_sample() {
  local sample=$1 code=$2
  local fq1=$FASTQ_DIR/Unknown_CM794-S01-001${code}_good_1.fq.gz
  local fq2=$FASTQ_DIR/Unknown_CM794-S01-001${code}_good_2.fq.gz
  local bam=$OUT/bam/${sample}.bam
  if [ -f "$bam" ]; then
    echo "[$(date '+%T')] skip (already present): $sample"
    return 0
  fi
  if [ ! -f "$fq1" ] || [ ! -f "$fq2" ]; then
    echo "[$(date '+%T')] [error] missing fastq: $sample ($code)" >&2
    return 1
  fi
  echo "[$(date '+%T')] alignment started: $sample ($code)"
  nice -n 10 bowtie2 -p $THREADS --very-sensitive -x "$GENESET" \
    -1 "$fq1" -2 "$fq2" 2> $LOG/step5_bt2_${sample}.log | \
    nice -n 10 samtools view -@ $SAM_THREADS -b - > "$bam"
  local rc=${PIPESTATUS[0]}
  if [ $rc -eq 0 ] && [ -s "$bam" ]; then
    echo "[$(date '+%T')] alignment finished: $sample"
  else
    echo "[$(date '+%T')] [error] alignment failed: $sample (bowtie2 rc=$rc)" >&2
    rm -f "$bam"
    return 1
  fi
}

echo "[$(date '+%F %T')] index ready, starting alignment (bowtie2 --very-sensitive, ${CONC} samples x ${THREADS} threads)"
# Sample mapping (cross-checked against archive size and read counts):
#   R0001=CK-1  R0002=CK-2  R0003=CK-3  R0004=HD-1  R0005=HD-2  R0006=HD-3
declare -a PAIRS=("CK-1 R0001" "CK-2 R0002" "CK-3 R0003" "HD-1 R0004" "HD-2 R0005" "HD-3 R0006")

FAILED=""
for pair in "${PAIRS[@]}"; do
  sample=${pair%% *}; code=${pair##* }
  map_one_sample "$sample" "$code" &
  # Cap concurrency: wait for a slot when the pool is full
  while [ "$(jobs -rp | wc -l)" -ge "$CONC" ]; do
    sleep 15
  done
done
wait
# Verify: all 6 BAM files must exist and be non-empty, otherwise abort
for pair in "${PAIRS[@]}"; do
  sample=${pair%% *}
  if [ ! -s "$OUT/bam/${sample}.bam" ]; then
    echo "[$(date '+%F %T')] [fatal] ${sample}.bam missing or empty" >&2
    FAILED="$FAILED $sample"
  fi
done
if [ -n "$FAILED" ]; then
  echo "[$(date '+%F %T')] alignment failed for these samples, aborting:$FAILED" >&2
  exit 1
fi

echo "[$(date '+%F %T')] alignment done, counting reads per gene (6 samples in parallel, 2 threads each)"
# Per-gene counts (samtools view + awk; light I/O task, run fully in parallel)
count_one_bam() {
  local bam=$1
  local sample=$(basename "$bam" .bam)
  nice -n 10 samtools view -@ 2 "$bam" | awk -v s="$sample" '
    and($2,0x100) {next}          # skip secondary
    and($2,0x800) {next}          # skip supplementary
    $3=="*" {next}                # unmapped
    {cnt[$3]++}
    END{for(g in cnt) print g"\t"cnt[g]}' > $OUT/${sample}.gene_counts.tsv
}
for bam in $OUT/bam/*.bam; do
  count_one_bam "$bam" &
done
wait

echo "[$(date '+%F %T')] Computing per-gene length and TPM"
# step5_tpm.py is located adaptively: 03_abundance first, then 01_pr_analysis
TPM_SCRIPT=$OUT/step5_tpm.py
[ -f "$TPM_SCRIPT" ] || TPM_SCRIPT=$PROJ/01_pr_analysis/step5_tpm.py
[ -f "$TPM_SCRIPT" ] || { echo "[fatal] step5_tpm.py not found" >&2; exit 1; }
python3 "$TPM_SCRIPT"

echo "[$(date '+%F %T')] Step5 complete"
