#!/bin/bash
# step4_tree.sh - PR phylogenetic tree (kelp OTUs + global references)
# env: phylo (trimal, iqtree) + Metagenomic tools on PATH first half
set -euo pipefail
BASE=$HOME/projects/kelpPR
A=$BASE/01_pr_analysis
G=$BASE/00_data/geneset
T=$A/tree
mkdir -p $T
cd $A

echo "=== [1/6] Dereplicate reference sequences ==="
cd-hit -i $G/uniprot_proteorhodopsins_all.fasta -o $T/pr_refs_derep.fasta -c 0.95 -aS 0.9 -d 0 -M 0 -T 8 -g 1 > $T/cdhit_pr.log 2>&1
cd-hit -i $G/uniprot_type1_rhodopsins_reviewed.fasta -o $T/t1_refs_derep.fasta -c 0.90 -aS 0.9 -d 0 -M 0 -T 8 -g 1 > $T/cdhit_t1.log 2>&1
echo "PR refs: $(grep -c '^>' $T/pr_refs_derep.fasta) | type-1 refs: $(grep -c '^>' $T/t1_refs_derep.fasta)"

echo "=== [2/6] Build tree input (rename labels) ==="
python3 step4_prep_tree_input.py $T/pr_refs_derep.fasta $T/t1_refs_derep.fasta $A/pr_otu_centroids_tree.fasta $T/tree_input.fasta

echo "=== [3/6] MAFFT alignment ==="
mafft --auto --thread 8 $T/tree_input.fasta > $T/tree_aln.fasta 2> $T/mafft_tree.log

echo "=== [4/6] trimAl trimming ==="
trimal -automated1 -in $T/tree_aln.fasta -out $T/tree_aln_trim.fasta

echo "=== [5/6] IQ-TREE（MFP + 1000 UFBoot）==="
cd $T
iqtree -s tree_aln_trim.fasta -m MFP -B 1000 -T 8 --prefix kelpPR_pr -seed 20260910

echo "=== [6/6] Results ==="
ls -la $T/kelpPR_pr.*
echo "=== ALL DONE (step4 tree) ==="
