#!/usr/bin/env python3
# =============================================================================
# Project : kelpPR (Dongshan Bay kelp aquaculture metagenome)
# Author  : Fei Zhang
# Date    : 2026-09-14
# Version : 1.0
# Input   : mag_pr_hits.tbl (hmmsearch tblout), mag_proteins/*.faa (lengths),
#           02_mag/gtdbtk/classify summaries, 03_pr_mag/pr_mag_mapping.tsv
# Output  : mag_pr_summary.tsv (MAG, n_PR, PR_proteins, mean_len, taxonomy,
#                              cross_validated_by_step7)
# Logic   : prodigal protein IDs are <contig>_<n>; MAG = protein id prefix up to
#           contig name... prodigal -p meta keeps original contig headers, so the
#           MAG of each hit is resolved via contig2bin.tsv from step7 output.
#           Length filter: 150-450 aa (parsed from .faa headers is not available;
#           use sequence lengths from all_mag_proteins.faa).
# =============================================================================
import csv, os
from collections import defaultdict

PROJ = os.environ.get("PROJECT_DIR", os.path.expanduser("~/projects/kelpPR"))
OUT  = os.path.join(PROJ, "05_mag_pr")

# sequence lengths
seqlen = {}
name, ln = None, 0
with open(os.path.join(OUT, "all_mag_proteins.faa")) as f:
    for line in f:
        if line.startswith(">"):
            if name: seqlen[name] = ln
            name, ln = line[1:].split()[0], 0
        else:
            ln += len(line.strip())
if name: seqlen[name] = ln

# contig -> bin
c2b = {}
with open(os.path.join(PROJ, "03_pr_mag/contig2bin.tsv")) as f:
    for line in f:
        c, b = line.rstrip("\n").split("\t")
        c2b[c] = b

# taxonomy
tax = {}
for sf in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"):
    p = os.path.join(PROJ, "02_mag/gtdbtk/classify", sf)
    if os.path.exists(p):
        with open(p) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                tax[row["user_genome"]] = row["classification"]

# step7 cross-validation: which MAGs got PR by read mapping
s7 = set()
p7 = os.path.join(PROJ, "03_pr_mag/pr_mag_mapping.tsv")
if os.path.exists(p7):
    with open(p7) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["MAG"] != "unbinned":
                s7.add(row["MAG"])

# parse hmmsearch hits: protein id = <contig>_<i>
hits = defaultdict(list)
with open(os.path.join(OUT, "mag_pr_hits.tbl")) as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.split()
        prot, evalue, score = p[0], float(p[4]), float(p[5])
        L = seqlen.get(prot, 0)
        if not (150 <= L <= 450): continue
        contig = prot.rsplit("_", 1)[0]
        mag = c2b.get(contig)
        if mag:
            hits[mag].append((prot, L, evalue))

rows = []
for mag in sorted(hits, key=lambda m: int(m.split(".")[1])):
    ps = hits[mag]
    rows.append((mag, len(ps), ",".join(x[0] for x in ps),
                 f"{sum(x[1] for x in ps)/len(ps):.0f}",
                 tax.get(mag, ""), "yes" if mag in s7 else "no"))

with open(os.path.join(OUT, "mag_pr_summary.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["MAG", "n_PR", "PR_proteins", "mean_len_aa", "GTDB_taxonomy", "step7_crossval"])
    w.writerows(rows)

print(f"MAGs with PR (direct evidence): {len(rows)}; of which cross-validated by step7: {sum(1 for r in rows if r[5]=='yes')}")
