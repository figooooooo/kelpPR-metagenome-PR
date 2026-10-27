#!/usr/bin/env python3
# =============================================================================
# Project : kelpPR (Dongshan Bay kelp aquaculture metagenome)
# Author  : Fei Zhang
# Date    : 2026-09-14
# Version : 1.0
# Input   : quant/bin_abundance_table.tab (metaWRAP quant_bins, salmon TPM)
#           02_mag/gtdbtk/classify/gtdbtk.{bac120,ar53}.summary.tsv
#           02_mag/checkm/qa.tsv
# Output  : mag_hd_vs_ck.tsv  (MAG, taxonomy, completeness, contamination,
#                              meanTPM_CK, meanTPM_HD, log2FC, exact MW p)
# Notes   : R0001-3 = CK, R0004-6 = HD. n=3v3 -> exact Mann-Whitney (20 combos),
#           min achievable p = 0.1; interpret descriptively.
# =============================================================================
import csv, os, re
from itertools import combinations

PROJ = os.environ.get("PROJECT_DIR", os.path.expanduser("~/projects/kelpPR"))
OUT  = os.path.join(PROJ, "04_mag_abundance")

def exact_mw_p(x, y):
    """Exact two-sided Mann-Whitney for n1=n2=3 via full enumeration."""
    vals = list(x) + list(y)
    n = len(x)
    obs_u = sum(1 for a in x for b in y if a > b) + 0.5 * sum(1 for a in x for b in y if a == b)
    idx = range(len(vals))
    cnt = tot = 0
    for comb in combinations(idx, n):
        us = sum(1 for i in comb for j in idx if j not in comb and vals[i] > vals[j]) \
             + 0.5 * sum(1 for i in comb for j in idx if j not in comb and vals[i] == vals[j])
        tot += 1
        if abs(us - n * n / 2) >= abs(obs_u - n * n / 2) - 1e-9:
            cnt += 1
    return cnt / tot

# taxonomy & quality
tax, qual = {}, {}
for sf in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"):
    p = os.path.join(PROJ, "02_mag/gtdbtk/classify", sf)
    if os.path.exists(p):
        with open(p) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                tax[row["user_genome"]] = row["classification"]
with open(os.path.join(PROJ, "02_mag/checkm/qa.tsv")) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        qual[row["Bin Id"]] = (row["Completeness"], row["Contamination"])

# abundance table
tab = os.path.join(OUT, "quant/bin_abundance_table.tab")
with open(tab) as f:
    header = f.readline().rstrip("\n").split("\t")
    # identify sample columns containing R000x
    col_group = {}
    for i, h in enumerate(header):
        m = re.search(r"R000([1-6])", h)
        if m:
            col_group[i] = "CK" if int(m.group(1)) <= 3 else "HD"
    bin_col = next(i for i, h in enumerate(header) if i not in col_group)
    rows = []
    for line in f:
        p = line.rstrip("\n").split("\t")
        name = p[bin_col]
        ck = [float(p[i]) for i in col_group if col_group[i] == "CK"]
        hd = [float(p[i]) for i in col_group if col_group[i] == "HD"]
        rows.append((name, ck, hd))

out_rows = []
for name, ck, hd in rows:
    mck, mhd = sum(ck) / 3, sum(hd) / 3
    l2fc = None
    if mck > 0 and mhd > 0:
        import math
        l2fc = math.log2(mhd / mck)
    p = exact_mw_p(hd, ck) if (mck > 0 or mhd > 0) else 1.0
    comp, cont = qual.get(name, ("", ""))
    out_rows.append((name, tax.get(name, ""), comp, cont,
                     f"{mck:.2f}", f"{mhd:.2f}",
                     f"{l2fc:+.2f}" if l2fc is not None else "inf" if mhd > 0 else "-inf" if mck > 0 else "0/0",
                     f"{p:.3f}"))

out_rows.sort(key=lambda r: -(float(r[4]) + float(r[5])))
with open(os.path.join(OUT, "mag_hd_vs_ck.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["MAG", "GTDB_taxonomy", "Completeness", "Contamination",
                "meanTPM_CK", "meanTPM_HD", "log2FC_HD_CK", "exact_MW_p"])
    w.writerows(out_rows)

hd_hi = sum(1 for r in out_rows if r[6] not in ("inf", "-inf", "0/0") and float(r[6]) >= 1)
ck_hi = sum(1 for r in out_rows if r[6] not in ("inf", "-inf", "0/0") and float(r[6]) <= -1)
print(f"MAGs: {len(out_rows)}; HD-enriched (log2FC>=1): {hd_hi}; CK-enriched (log2FC<=-1): {ck_hi}")
