#!/bin/bash
# step2_pipeline.sh - PR filtering, functional-site annotation, OTU clustering
# env: Metagenomic (mafft, cd-hit, python3)
set -euo pipefail
BASE=$HOME/projects/kelpPR
A=$BASE/01_pr_analysis
G=$BASE/00_data/geneset
cd $A

echo "=== [1/5] Extract green/blue-annotated reference sequences ==="
python3 - "$G/uniprot_proteorhodopsins_all.fasta" "$A/calib_refs.fasta" <<'EOF'
import sys
def read_fasta(path):
    order, seqs, name, buf = [], {}, None, []
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('>'):
            if name is not None: seqs[name] = ''.join(buf)
            name = line[1:]; order.append(name); buf = []
        else:
            buf.append(line.strip())
    if name is not None: seqs[name] = ''.join(buf)
    return order, seqs
order, seqs = read_fasta(sys.argv[1])
n = 0
with open(sys.argv[2], 'w') as f:
    for h in order:
        if 'green' in h.lower() or 'blue' in h.lower():
            f.write('>%s\n' % h)
            s = seqs[h]
            for i in range(0, len(s), 60): f.write(s[i:i+60] + '\n')
            n += 1
print('calibration refs: %d' % n)
EOF

echo "=== [2/5] MAFFT alignment (candidates + calibration references) ==="
cat pr_candidates_len.fasta calib_refs.fasta > calib_input.fasta
mafft --auto --thread 24 calib_input.fasta > calib_aln.fasta 2> mafft_calib.log
echo "alignment done: $(grep -c '^>' calib_aln.fasta) seqs"

echo "=== [3/5] Filtering + functional-site annotation ==="
python3 step2_pr_filter.py calib_aln.fasta hmmsearch_pf01036.tbl $A

echo "=== [4/5] CD-HIT 95% OTU clustering (parameters from Li 2024) ==="
cd-hit -i pr_filtered.fasta -o pr_otus_95.fasta -c 0.95 -aS 0.9 -d 0 -M 0 -T 8 -g 1 > cdhit_otu.log 2>&1

echo "=== [5/5] OTU abundance table ==="
python3 step3_otu_table.py pr_otus_95.fasta.clstr pr_filtered.fasta pr_annotation_all.tsv $A

echo "=== ALL DONE (step2) ==="
