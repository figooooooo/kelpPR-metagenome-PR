#!/usr/bin/env python3
# =============================================================================
# Project : kelpPR (Dongshan Bay kelp aquaculture metagenome)
# Author  : Fei Zhang
# Date    : 2026-09-14
# Version : 1.0
# Input   : pr_vs_coassembly.blast6, contig2bin.tsv,
#           02_mag/gtdbtk/classify/gtdbtk.{bac120,ar53}.summary.tsv
#           01_pr_analysis/pr_otus_95.fasta.clstr (PR gene -> 95% OTU centroid)
# Output  : pr_mag_mapping.tsv
#   PR_gene, PR_OTU, best_hit_contig, pident, qcovs, MAG, GTDB_taxonomy
# Logic   : best blast hit per PR (highest bitscore); require pident>=95 & qcovs>=80;
#           contig in MAG -> MAG + taxonomy; else "unbinned".
# =============================================================================
import csv, os, re

PROJ = os.environ.get("PROJECT_DIR", os.path.expanduser("~/projects/kelpPR"))
OUT  = os.path.join(PROJ, "03_pr_mag")

# contig -> bin
c2b = {}
with open(os.path.join(OUT, "contig2bin.tsv")) as f:
    for line in f:
        ctg, b = line.rstrip("\n").split("\t")
        c2b[ctg] = b

# bin -> taxonomy
tax = {}
for sf in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"):
    p = os.path.join(PROJ, "02_mag/gtdbtk/classify", sf)
    if not os.path.exists(p):
        continue
    with open(p) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            tax[row["user_genome"]] = row["classification"]

# PR gene -> OTU centroid (from cd-hit .clstr)
gene2otu = {}
clstr = os.path.join(PROJ, "01_pr_analysis/pr_otus_95.fasta.clstr")
cur = None
if os.path.exists(clstr):
    with open(clstr) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                cur = line.split()[1]
            elif line:
                m = re.search(r">(\S+?)\.\.\.", line)
                if m:
                    g = m.group(1)
                    if line.endswith("*"):
                        gene2otu.setdefault(g, g)  # centroid placeholder
                    cur_id = None
    # second pass: map members to centroid id
    cent = None
    with open(clstr) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                cent = None
            elif line:
                m = re.search(r">(\S+?)\.\.\.", line)
                if m and line.endswith("*"):
                    cent = m.group(1)
                    gene2otu[cent] = cent
                elif m and cent:
                    gene2otu[m.group(1)] = cent

# best hit per PR
best = {}
with open(os.path.join(OUT, "pr_vs_coassembly.blast6")) as f:
    for line in f:
        q, s, pident, length, qlen, slen, qcovs, ev, bits = line.rstrip("\n").split("\t")
        pident, qcovs, bits = float(pident), float(qcovs), float(bits)
        if pident < 95 or qcovs < 80:
            continue
        if q not in best or bits > best[q][2]:
            best[q] = (s, pident, bits, qcovs)

rows = []
for q in sorted(best):
    s, pident, bits, qcovs = best[q]
    mag = c2b.get(s, "unbinned")
    t = tax.get(mag, "") if mag != "unbinned" else ""
    rows.append((q, gene2otu.get(q, ""), s, f"{pident:.1f}", f"{qcovs:.0f}", mag, t))

with open(os.path.join(OUT, "pr_mag_mapping.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["PR_gene", "PR_OTU", "coasm_contig", "pident", "qcovs", "MAG", "GTDB_taxonomy"])
    w.writerows(rows)

n_mapped = sum(1 for r in rows if r[5] != "unbinned")
print(f"PR with valid hit: {len(rows)}; mapped to MAG: {n_mapped}; unbinned: {len(rows)-n_mapped}")
