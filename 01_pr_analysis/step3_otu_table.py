#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step3_otu_table.py - OTU (95% aa identity, cd-hit) abundance table + tree-ready centroid fasta
Usage: python3 step3_otu_table.py <clstr> <filtered_fasta> <annotation_tsv> <outdir>
"""
import sys
from collections import Counter

CLSTR, FFASTA, TSVO, OUTDIR = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

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

anno = {}
for i, line in enumerate(open(TSVO)):
    if i == 0: continue
    f = line.rstrip('\n').split('\t')
    if len(f) >= 11:
        anno[f[0]] = {'sample': f[1], 'group': f[2], 'spectral': f[7], 'motif': f[8]}

forder, fseqs = read_fasta(FFASTA)
bygene = {h.split()[0]: h for h in forder}

# parse cd-hit .clstr
otus, cur, cent = [], None, None
for line in open(CLSTR):
    line = line.rstrip('\n')
    if line.startswith('>Cluster'):
        if cur is not None: otus.append((cent, cur))
        cur, cent = [], None
    else:
        # format: "0\t238aa, >geneset000234664 protein=... at 1:238:1:238 *"
        gid = line.split('>', 1)[1].split()[0].rstrip('.')
        cur.append(gid)
        if line.rstrip().endswith('*'): cent = gid
if cur is not None: otus.append((cent, cur))

samples = sorted({a['sample'] for a in anno.values()} |
                 {'CK-1', 'CK-2', 'CK-3', 'HD-1', 'HD-2', 'HD-3'})
with open(OUTDIR + '/pr_otu_table.tsv', 'w') as f:
    f.write('OTU_id\tcentroid\tsize\t' + '\t'.join(samples) + '\tconsensus_spectral\tconsensus_motif\n')
    for i, (cent, members) in enumerate(sorted(otus, key=lambda x: -len(x[1])), 1):
        counts = Counter(anno[g]['sample'] for g in members if g in anno)
        sp = Counter(anno[g]['spectral'] for g in members if g in anno).most_common(1)
        mo = Counter(anno[g]['motif'] for g in members if g in anno).most_common(1)
        f.write('OTU%03d\t%s\t%d\t%s\t%s\t%s\n' % (
            i, cent, len(members),
            '\t'.join(str(counts.get(s, 0)) for s in samples),
            sp[0][0] if sp else 'NA', mo[0][0] if mo else 'NA'))

# tree-ready centroids with informative labels
with open(OUTDIR + '/pr_otu_centroids_tree.fasta', 'w') as f:
    for i, (cent, members) in enumerate(sorted(otus, key=lambda x: -len(x[1])), 1):
        h = bygene.get(cent, '')
        if not h: continue
        seq = fseqs[h]
        a = anno.get(cent, {})
        tag = {'green': 'GPR', 'blue': 'BPR', 'unknown': 'UK'}.get(a.get('spectral', 'UK'), 'UK')
        groups = set(anno[g]['group'] for g in members if g in anno)
        label = 'KelpOTU%03d|%s|%s|n%d' % (i, tag, '-'.join(sorted(groups)) or '?', len(members))
        f.write('>%s\n' % label)
        for j in range(0, len(seq), 60): f.write(seq[j:j+60] + '\n')

print('[done] %d OTUs written (table + centroids)' % len(otus))
for i, (cent, members) in enumerate(sorted(otus, key=lambda x: -len(x[1]))[:10], 1):
    print('  OTU%03d size=%d centroid=%s' % (i, len(members), cent))
