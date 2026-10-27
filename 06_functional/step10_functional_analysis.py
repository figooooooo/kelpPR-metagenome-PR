#!/usr/bin/env python3
# =============================================================================
# Project : kelpPR (Dongshan Bay kelp aquaculture metagenome)
# Author  : Fei Zhang
# Date    : 2026-09-14
# Version : 1.0
# Input   : 00_data/annotation/kegg_anno/kegg.kegg_pathway_level{2,3}.relabundance.xls
#           00_data/annotation/kegg_anno/kegg.KO.gene_list.xls
#           00_data/annotation/CAZy_anno/CAZy.{Class,Family}.relabundance.xls
#           00_data/annotation/CAZy_anno/CAZy.anno_result.xls
#           00_data/annotation/Nr_anno/nr.taxonomy.abundance.xls
#           03_abundance/pr_tpm.tsv  (139 PR genes x 6 samples)
# Output  : 06_functional/kegg_l2_diff.tsv, kegg_l3_diff.tsv,
#           cazy_class_diff.tsv, cazy_family_diff.tsv,
#           kelp_polysaccharide_families.tsv, pr_ko_annotation.tsv,
#           pr_cazy_coupling.tsv, nr_genus_top.tsv
# Notes   : CK-1..3 = control, HD-1..3 = kelp farm. n=3v3 -> exact
#           Mann-Whitney (20 combos, min p = 0.1), identical implementation to
#           step8_compare.py; BH-FDR across each tested family. Relabundance
#           tables are proportions of annotated reads (Biomarker delivery).
# =============================================================================
import os
from itertools import combinations

PROJ = os.environ.get("PROJECT_DIR", os.path.expanduser("~/projects/kelpPR"))
ANN = os.path.join(PROJ, "00_data/annotation")
OUT = os.path.join(PROJ, "06_functional")
os.makedirs(OUT, exist_ok=True)

CK = ["CK-1", "CK-2", "CK-3"]
HD = ["HD-1", "HD-2", "HD-3"]
SAMPLES = CK + HD


def exact_mw_p(x, y):
    """Exact two-sided Mann-Whitney for n1=n2=3 via full enumeration (as step8)."""
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


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR. Returns adjusted q-values in original order."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        val = min(prev, pvals[i] * m / (rank + 1))
        q[i] = val
        prev = val
    return q


def read_table(path):
    rows = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                rows.append(p)
    return header, rows


def diff_test(in_path, out_path, name_col, extra_cols, ck_mean_out="mean_CK", hd_mean_out="mean_HD"):
    """Generic relabundance differential test -> tsv."""
    header, rows = read_table(in_path)
    si = {s: header.index(s) for s in SAMPLES}
    results = []
    for p in rows:
        try:
            ckv = [float(p[si[s]]) for s in CK]
            hdv = [float(p[si[s]]) for s in HD]
        except (ValueError, IndexError):
            continue
        mck, mhd = sum(ckv) / 3, sum(hdv) / 3
        ratio = (mhd + 1e-12) / (mck + 1e-12)
        pv = exact_mw_p(ckv, hdv)
        results.append({
            "name": p[0], "extras": [p[i] if i < len(p) else "" for i in extra_cols],
            "mck": mck, "mhd": mhd, "ratio": ratio,
            "l2fc": ratio if ratio == 0 else (0 if mck == 0 and mhd == 0 else ratio),
            "p": pv, "dir": "HD" if mhd > mck else "CK",
        })
    qs = bh_fdr([r["p"] for r in results])
    with open(out_path, "w") as fo:
        fo.write("\t".join([name_col] + [header[i] for i in extra_cols if i < len(header)]
                           + [ck_mean_out, hd_mean_out, "ratio_HD_CK", "log2FC",
                              "direction", "p_exactMW", "q_BH"]) + "\n")
        for r, q in zip(results, qs):
            import math
            l2fc = math.log2(r["ratio"]) if r["mck"] > 0 and r["mhd"] > 0 else ("inf_HD" if r["mhd"] > 0 else ("neg_inf" if r["mck"] > 0 else "NA"))
            fo.write("\t".join([r["name"]] + r["extras"]
                               + ["%.6g" % r["mck"], "%.6g" % r["mhd"], "%.4f" % r["ratio"],
                                  str(l2fc), r["dir"], "%.4f" % r["p"], "%.4f" % q]) + "\n")
    n_sig = sum(1 for r, q in zip(results, qs) if r["p"] < 0.05)
    print("[step10] %s: %d features, p<0.05 unadjusted: %d (min p = 0.1 by design)"
          % (os.path.basename(out_path), len(results), n_sig))
    return results


# ---------- 1. KEGG L2 / L3 differential ----------
diff_test(os.path.join(ANN, "kegg_anno/kegg.kegg_pathway_level2.relabundance.xls"),
          os.path.join(OUT, "kegg_l2_diff.tsv"), "pathway", [1], )
l3 = diff_test(os.path.join(ANN, "kegg_anno/kegg.kegg_pathway_level3.relabundance.xls"),
               os.path.join(OUT, "kegg_l3_diff.tsv"), "pathway", [1, 2, 3])

# ---------- 2. CAZy Class / Family differential ----------
diff_test(os.path.join(ANN, "CAZy_anno/CAZy.Class.relabundance.xls"),
          os.path.join(OUT, "cazy_class_diff.tsv"), "class", [1])
fam = diff_test(os.path.join(ANN, "CAZy_anno/CAZy.Family.relabundance.xls"),
                os.path.join(OUT, "cazy_family_diff.tsv"), "family", [1])

# kelp-polysaccharide focus families (laminarin GH16/GH17, alginate PL6/PL7/PL17, fucoidan GH29/GH95/GH107/GH168..., agar GH16/GH86/GH118, CBMs)
FOCUS = ["GH16", "GH17", "GH3", "GH5", "GH30", "GH29", "GH95", "GH107", "GH117", "GH168", "GH169",
         "PL1", "PL6", "PL7", "PL14", "PL17", "PL18", "CBM16", "CBM32", "CBM35", "CBM48", "AA7"]
with open(os.path.join(OUT, "kelp_polysaccharide_families.tsv"), "w") as fo:
    fo.write("family\tactivity\tmean_CK\tmean_HD\tratio_HD_CK\tp_exactMW\tdirection\n")
    for r in fam:
        if r["name"] in FOCUS:
            fo.write("%s\t%s\t%.6g\t%.6g\t%.4f\t%.4f\t%s\n"
                     % (r["name"], r["extras"][0], r["mck"], r["mhd"], r["ratio"], r["p"], r["dir"]))
print("[step10] kelp_polysaccharide_families.tsv written (%d focus families)" % len(FOCUS))

# ---------- 3. PR gene KO annotation (cross-check) ----------
pr_ids, pr_tpm = [], {}
with open(os.path.join(PROJ, "03_abundance/pr_tpm.tsv")) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        pr_ids.append(p[0])
        pr_tpm[p[0]] = [float(v) for v in p[1:7]]
pr_set = set(pr_ids)
ko_map = {}
with open(os.path.join(ANN, "kegg_anno/kegg.KO.gene_list.xls"), encoding="utf-8", errors="replace") as fh:
    fh.readline()
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 3:
            continue
        for g in p[2].split(";"):
            if g in pr_set:
                ko_map.setdefault(g, []).append((p[0], p[1]))
with open(os.path.join(OUT, "pr_ko_annotation.tsv"), "w") as fo:
    fo.write("pr_gene\tKO_id\tKO_function\tTPM_CK_mean\tTPM_HD_mean\n")
    for g in pr_ids:
        for ko, kfn in ko_map.get(g, [("NA", "no KO assigned")]):
            ck_m = sum(pr_tpm[g][:3]) / 3
            hd_m = sum(pr_tpm[g][3:]) / 3
            fo.write("%s\t%s\t%s\t%.3f\t%.3f\n" % (g, ko, kfn, ck_m, hd_m))
n_ko = len([g for g in pr_ids if g in ko_map])
print("[step10] PR genes with KO: %d / %d -> pr_ko_annotation.tsv" % (n_ko, len(pr_ids)))

# ---------- 4. PR x CAZy sample-level coupling (Spearman, n=6) ----------
def rank(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r

def spearman(x, y):
    rx, ry = rank(x), rank(y)
    mx, my = sum(rx) / 6, sum(ry) / 6
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else 0.0

pr_sum = [sum(pr_tpm[g][i] for g in pr_ids) for i in range(6)]
cazy_fam_hdr, cazy_fam_rows = read_table(os.path.join(ANN, "CAZy_anno/CAZy.Family.relabundance.xls"))
si = {s: cazy_fam_hdr.index(s) for s in SAMPLES}
coupling = []
for p in cazy_fam_rows:
    if not p[0]:
        continue
    vals = [float(p[si[s]]) for s in SAMPLES]
    coupling.append((p[0], spearman(pr_sum, vals), p[1]))
with open(os.path.join(OUT, "pr_cazy_coupling.tsv"), "w") as fo:
    fo.write("family\tspearman_prTpm_vs_family\tactivity\n")
    for name, rho, act in sorted(coupling, key=lambda t: -abs(t[1])):
        fo.write("%s\t%.3f\t%s\n" % (name, rho, act[:80]))
focus_rho = {n: r for n, r, _ in coupling if n in FOCUS}
print("[step10] PR-CAZy coupling (Spearman, n=6): " +
      ", ".join("%s=%.2f" % (n, focus_rho[n]) for n in FOCUS if n in focus_rho))

# ---------- 5. NR genus-level top table ----------
gn_hdr, gn_rows = read_table(os.path.join(ANN, "Nr_anno/nr.taxonomy.abundance.xls"))
gsi = {s: gn_hdr.index(s) for s in SAMPLES}
genera = []
for p in gn_rows:
    parts = [x.strip() for x in p[0].split(";")]
    if len(parts) >= 6 and parts[5].startswith("g__") and parts[5] != "g__Unclassified":
        try:
            vals = [float(p[gsi[s]]) for s in SAMPLES]
        except (ValueError, IndexError):
            continue
        genera.append((parts[5][3:], parts[1][3:] if len(parts) > 1 else "",
                        vals))
genera.sort(key=lambda t: -sum(t[2]) / 6)
with open(os.path.join(OUT, "nr_genus_top.tsv"), "w") as fo:
    fo.write("genus\tphylum\t" + "\t".join(SAMPLES) + "\tmean_CK\tmean_HD\tp_exactMW\n")
    for gname, ph, vals in genera[:40]:
        pv = exact_mw_p(vals[:3], vals[3:])
        fo.write("%s\t%s\t" % (gname, ph) + "\t".join("%.1f" % v for v in vals)
                 + "\t%.1f\t%.1f\t%.4f\n" % (sum(vals[:3]) / 3, sum(vals[3:]) / 3, pv))
print("[step10] nr_genus_top.tsv written (top 40 of %d genera)" % len(genera))

print("[step10] ALL DONE -> " + OUT)
