#!/usr/bin/env python3
"""Step5 post-processing: TPM calculation + PR abundance GPM normalisation + CK vs HD group test
Inputs (all under ~/projects/kelpPR/):
  00_data/geneset/geneset_gene.fasta     nucleotide gene catalogue (headers carry len=)
  00_data/geneset/scg_ko_map.tsv         gene ID<TAB>KO<TAB>name (recA/rpoB/rplB)
  01_pr_analysis/pr_pass_ids.txt         139 PASS PR gene IDs
  03_abundance/<sample>.gene_counts.tsv  per-gene read counts (written by the step5 driver)
Outputs (03_abundance/):
  gene_tpm_all.tsv          TPM matrix for all genes (samples x genes)
  pr_tpm.tsv                TPM of the 139 PR genes
  pr_gpm_summary.tsv        per-sample PR-GPM, PR-TPM share, median SCG-TPM
  pr_group_test.tsv         CK vs HD exact Mann-Whitney test (3 vs 3, full enumeration)
"""
import os, sys, math
from collections import defaultdict

PROJ = os.environ.get("PROJECT_DIR", os.path.expanduser("~/projects/kelpPR"))
GA = os.path.join(PROJ, "00_data/geneset/geneset_gene.fasta")
SCG = os.path.join(PROJ, "00_data/geneset/scg_ko_map.tsv")
PR_IDS = os.path.join(PROJ, "01_pr_analysis/pr_pass_ids.txt")
OUT = os.path.join(PROJ, "03_abundance")

# ---- sample grouping (delivery file order -> sample name, see 00_data/fastq mapping) ----
GROUP = {"CK-1": "CK", "CK-2": "CK", "CK-3": "CK", "HD-1": "HD", "HD-2": "HD", "HD-3": "HD"}

# ---- gene lengths ----
glen = {}
with open(GA) as fh:
    for line in fh:
        if line.startswith(">"):
            parts = line[1:].split()
            gid = parts[0]
            for p in parts[1:]:
                if p.startswith("len="):
                    glen[gid] = int(p[4:])
                    break

# ---- per-sample gene counts ----
samples = []
counts = {}  # sample -> {gene: count}
for f in sorted(os.listdir(OUT)):
    if f.endswith(".gene_counts.tsv"):
        s = f.replace(".gene_counts.tsv", "")
        samples.append(s)
        counts[s] = {}
        with open(os.path.join(OUT, f)) as fh:
            for line in fh:
                g, c = line.rstrip("\n").split("\t")
                counts[s][g] = int(c)
samples = [s for s in ["CK-1","CK-2","CK-3","HD-1","HD-2","HD-3"] if s in samples] or samples

# ---- TPM ----
tpm = {s: {} for s in samples}
for s in samples:
    rpk_sum = 0.0
    rpk = {}
    for g, c in counts[s].items():
        L = glen.get(g, 1)
        r = c / (L / 1000.0)
        rpk[g] = r
        rpk_sum += r
    for g, r in rpk.items():
        tpm[s][g] = r / rpk_sum * 1e6 if rpk_sum > 0 else 0.0

# ---- TPM matrix for all genes (only genes with counts, to save space) ----
all_genes = sorted(set().union(*[set(counts[s]) for s in samples]))
with open(os.path.join(OUT, "gene_tpm_all.tsv"), "w") as fh:
    fh.write("gene_id\t" + "\t".join(samples) + "\n")
    for g in all_genes:
        fh.write(g + "\t" + "\t".join(f"{tpm[s].get(g,0):.4f}" for s in samples) + "\n")

# ---- PR TPM ----
pr_ids = [l.strip() for l in open(PR_IDS) if l.strip()]
pr_set = set(pr_ids)
with open(os.path.join(OUT, "pr_tpm.tsv"), "w") as fh:
    fh.write("gene_id\t" + "\t".join(samples) + "\tsum\n")
    for g in pr_ids:
        if any(g in counts[s] for s in samples):
            fh.write(g + "\t" + "\t".join(f"{tpm[s].get(g,0):.4f}" for s in samples) +
                     f"\t{sum(tpm[s].get(g,0) for s in samples):.4f}\n")

# ---- SCG TPM ----
scg_ko = {}
for line in open(SCG):
    g, ko, name = line.rstrip("\n").split("\t")
    scg_ko[g] = ko

# ---- GPM: PR sum-RPK / SCG mean-RPK (RPK is implicit in the TPM normalisation, so the TPM ratio is equivalent) ----
rows = []
for s in samples:
    pr_sum = sum(tpm[s].get(g, 0) for g in pr_ids)
    ko_sums = defaultdict(list)
    for g, ko in scg_ko.items():
        if g in tpm[s] and tpm[s][g] > 0:
            ko_sums[ko].append(tpm[s][g])
    scg_mean = {ko: sum(v)/len(v) for ko, v in ko_sums.items()}
    scg_med = sorted(scg_mean.values())[len(scg_mean)//2] if scg_mean else 0
    gpm = pr_sum / scg_med if scg_med > 0 else float("nan")
    rows.append((s, GROUP.get(s, "?"), pr_sum, gpm, scg_mean.get("K03553", 0), scg_mean.get("K03043", 0), scg_mean.get("K02886", 0)))

with open(os.path.join(OUT, "pr_gpm_summary.tsv"), "w") as fh:
    fh.write("sample\tgroup\tPR_TPM_sum\tPR_GPM\trecA_meanTPM\trpoB_meanTPM\trplB_meanTPM\n")
    for r in rows:
        fh.write("\t".join([r[0], r[1], f"{r[2]:.2f}", f"{r[3]:.4f}", f"{r[4]:.2f}", f"{r[5]:.2f}", f"{r[6]:.2f}"]) + "\n")

# ---- group test (exact Mann-Whitney for 3 vs 3 = 20 permutations) ----
def exact_mw(a, b):
    """Return (U, one-sided p in the requested direction, two-sided p)."""
    n, m = len(a), len(b)
    allv = [(v, 0) for v in a] + [(v, 1) for v in b]
    allv.sort()
    # ranks (no ties assumed)
    ranks = {}
    for i, (v, _) in enumerate(allv):
        ranks.setdefault(v, []).append(i + 1)
    r = {tuple(): 0}
    U1 = sum(sum(ranks[v][0] for _ in [0]) for v in a) if False else None
    # simplified: compute U1 = R1 - n(n+1)/2 directly
    R1 = 0
    i = 0
    for v, g in allv:
        i += 1
        if g == 0:
            R1 += i
    U1 = R1 - n*(n+1)//2
    # exact p by full enumeration (C(6,3) = 20)
    from itertools import combinations
    combined = a + b
    obs = sum(a)
    cnt = 0; total = 0; ge = 0; le = 0
    for comb in combinations(range(n+m), n):
        ssum = sum(combined[i] for i in comb)
        total += 1
        if ssum >= obs: ge += 1
        if ssum <= obs: le += 1
    p_two = min(1.0, 2 * min(ge, le) / total)
    return U1, min(ge,le)/total, p_two, ge, le, total

for metric, idx in [("PR_TPM_sum", 2), ("PR_GPM", 3)]:
    ck = [r[idx] for r in rows if r[1] == "CK"]
    hd = [r[idx] for r in rows if r[1] == "HD"]
    if ck and hd:
        U1, p1, p2, ge, le, tot = exact_mw(ck, hd)
        with open(os.path.join(OUT, "pr_group_test.tsv"), "a") as fh:
            fh.write(f"{metric}\tCK={ck}\tHD={hd}\tU={U1}\tp_greater={ge}/{tot}\tp_less={le}/{tot}\tp_two={p2:.4f}\n")

print("step5_tpm.py done: wrote gene_tpm_all.tsv / pr_tpm.tsv / pr_gpm_summary.tsv / pr_group_test.tsv")
