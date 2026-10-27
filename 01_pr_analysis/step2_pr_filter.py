#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step2_pr_filter.py (v2) - PR candidate filtering and functional-site annotation
Changes vs v1:
  - fixed consensus() bug (now indexes alignment sequences)
  - TM check replaced by hmmsearch domain coverage (>=80% of 224-aa Bac_rhodopsin model),
    since Kyte-Doolittle hydropathy is unreliable for charged PR helices
  - Phobius/DeepTMHMM web verification of 7-TM = pre-submission QC step
Inputs:
  calib_aln.fasta          - MAFFT alignment (candidates + green/blue UniProt PRs)
  hmmsearch_pf01036.tbl    - hmmsearch tblout (E-values)
  hmmsearch_pf01036.domtbl - hmmsearch domtblout (model coverage)
Outputs:
  pr_filtered.fasta, pr_annotation_all.tsv
"""
import sys, re
from collections import Counter, defaultdict

CAND_ALN = sys.argv[1] if len(sys.argv) > 1 else 'calib_aln.fasta'
TBL = sys.argv[2] if len(sys.argv) > 2 else 'hmmsearch_pf01036.tbl'
DOMTBL = sys.argv[3] if len(sys.argv) > 3 else 'hmmsearch_pf01036.domtbl'
OUTDIR = sys.argv[4] if len(sys.argv) > 4 else '.'

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

def ungapped_index(seq, col):
    return sum(1 for ch in seq[:col] if ch != '-')

def col_of_ungapped(seq, idx):
    n = 0
    for c, ch in enumerate(seq):
        if ch != '-':
            if n == idx: return c
            n += 1
    return None

def consensus(seqlist, col):
    c = Counter(s[col] for s in seqlist if col < len(s) and s[col] != '-')
    return c.most_common(1)[0][0] if c else '-'

def residue(seq, col):
    return seq[col] if col is not None and col < len(seq) and seq[col] != '-' else '-'

# ---------- main ----------
order, aln = read_fasta(CAND_ALN)
L = max(len(s) for s in aln.values())

cands, green, blue = [], [], []
for h in order:
    hl = h.lower()
    if 'green' in hl: green.append(h)
    elif 'blue' in hl: blue.append(h)
    else: cands.append(h)

print('[info] alignment: %d cols | cands=%d green=%d blue=%d' % (L, len(cands), len(green), len(blue)))

# --- E-values & model coverage ---
evals, covs = {}, {}
for line in open(TBL):
    if line.startswith('#'): continue
    f = line.split()
    if len(f) > 4: evals.setdefault(f[0], f[4])
for line in open(DOMTBL):
    if line.startswith('#'): continue
    f = line.split()
    if len(f) > 14:
        gid, qlen, hfrom, hto = f[0], int(f[5]), int(f[15]), int(f[16])
        covs[gid] = (hto - hfrom + 1) / qlen

# --- 1. retinal-binding Lys column (~231) ---
pr_refs = [aln[h] for h in green + blue]
kcols = []
for c in range(L):
    if c < 0.5 * L: continue
    ks = sum(1 for s in pr_refs if c < len(s) and s[c] == 'K')
    if ks / max(len(pr_refs), 1) >= 0.80:
        kcols.append(c)
kcol = kcols[-1] if kcols else None
print('[K231] candidate K columns: %s -> using %s (%.0f%% of alignment)' % (kcols, kcol, 100.0*kcol/L if kcol is not None else -1))

# --- 2. spectral position 105 (green M/L vs blue Q) ---
gseqs = [aln[h] for h in green]
bseqs = [aln[h] for h in blue]
c105_candidates = []
for c in range(L):
    gcons, bcons = consensus(gseqs, c), consensus(bseqs, c)
    if gcons in 'ML' and bcons == 'Q':
        c105_candidates.append(c)
if c105_candidates:
    c105 = min(c105_candidates, key=lambda c: abs(c - 0.42 * L))
else:
    c105 = None
print('[pos105] matching columns: %s -> using %s (expect ~0.42L=%.0f)' % (c105_candidates, c105, 0.42 * L))

# --- 3. positions 97/101/108 (-8/-4/+3 relative to 105) ---
c97 = c101 = c108 = None
if c105 is not None:
    votes97, votes101, votes108 = [], [], []
    for h in green:
        s = aln[h]
        pi = ungapped_index(s, c105)
        if pi is None or pi < 8: continue
        a, b, d = col_of_ungapped(s, pi-8), col_of_ungapped(s, pi-4), col_of_ungapped(s, pi+3)
        if a is not None: votes97.append(a)
        if b is not None: votes101.append(b)
        if d is not None: votes108.append(d)
    def modal(v): return Counter(v).most_common(1)[0][0] if v else None
    c97, c101, c108 = modal(votes97), modal(votes101), modal(votes108)
print('[97/101/108] columns: %s/%s/%s' % (c97, c101, c108))

gmotifs = Counter(residue(aln[h], c97) + residue(aln[h], c101) + residue(aln[h], c108) for h in green)
bmotifs = Counter(residue(aln[h], c97) + residue(aln[h], c101) + residue(aln[h], c108) for h in blue)
print('[verify] green refs 97/101/108 motifs: %s' % dict(gmotifs))
print('[verify] blue  refs 97/101/108 motifs: %s' % dict(bmotifs))

# --- full sequences ---
forder, fseqs = read_fasta('../01_pr_analysis/pr_candidates_len.fasta')

# --- annotate candidates ---
rows, pass_list = [], []
for h in cands:
    gene_id = h.split()[0]
    m = re.search(r'protein=(CK|HD)-(\d)', h)
    group = m.group(1) if m else '?'
    sample = '%s-%s' % (m.group(1), m.group(2)) if m else '?'
    full = fseqs.get(h, fseqs.get(gene_id, ''))
    resk = residue(aln[h], kcol)
    r105 = residue(aln[h], c105)
    spectral = 'green' if r105 in 'ML' else ('blue' if r105 == 'Q' else 'unknown')
    motif = residue(aln[h], c97) + residue(aln[h], c101) + residue(aln[h], c108)
    cov = covs.get(gene_id, 0.0)
    ok = (resk == 'K') and (200 <= len(full) <= 350) and (cov >= 0.80)
    rows.append((gene_id, sample, group, len(full), evals.get(gene_id, 'NA'),
                 resk, r105, spectral, motif, round(cov, 3), 'PASS' if ok else 'FAIL'))
    if ok: pass_list.append((h, full))

with open(OUTDIR + '/pr_annotation_all.tsv', 'w') as f:
    f.write('gene_id\tsample\tgroup\tlen\tEvalue\tK231\tres105\tspectral\tmotif_97_101_108\tmodel_cov\tfilter\n')
    for r in rows: f.write('\t'.join(str(x) for x in r) + '\n')
with open(OUTDIR + '/pr_filtered.fasta', 'w') as f:
    for h, seq in pass_list:
        f.write('>%s\n' % h.strip())
        for i in range(0, len(seq), 60): f.write(seq[i:i+60] + '\n')

# --- summary ---
n = len(rows)
npass = sum(1 for r in rows if r[10] == 'PASS')
by = defaultdict(Counter)
for r in rows:
    if r[10] == 'PASS': by[r[2]][r[7]] += 1
print('\n===== SUMMARY =====')
print('candidates: %d | PASS: %d | FAIL: %d' % (n, npass, n - npass))
for g in sorted(by):
    print('  %s: %s (total %d)' % (g, dict(by[g]), sum(by[g].values())))
sp = Counter(r[7] for r in rows if r[10] == 'PASS')
mo = Counter(r[8] for r in rows if r[10] == 'PASS')
covd = Counter('>=0.9' if r[9] >= 0.9 else ('0.8-0.9' if r[9] >= 0.8 else '<0.8') for r in rows)
print('spectral types (PASS): %s' % dict(sp))
print('motifs (PASS): %s' % dict(mo))
print('model coverage dist (all): %s' % dict(covd))
kfail = sum(1 for r in rows if r[5] != 'K')
print('K231 != K: %d | cov<0.8: %d | len out: %d' % (
    kfail, sum(1 for r in rows if r[9] < 0.8), sum(1 for r in rows if not 200 <= r[3] <= 350)))
