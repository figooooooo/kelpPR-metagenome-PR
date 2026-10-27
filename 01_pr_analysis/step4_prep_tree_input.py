#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step4_prep_tree_input.py - build tree reference set with clean labels
- UniProt PRs: green-labeled (GPR), blue-labeled (BPR), others (PR, derep by user beforehand)
- Reviewed type-1 rhodopsins: functional class labels (bR/xR/HR/sR/ChR/ActR/ArchR/...)
- combined with kelp OTU centroids
Usage: python3 step4_prep_tree_input.py <pr_all_fasta> <type1_fasta> <otu_centroid_fasta> <outfile>
"""
import sys, re
from collections import Counter

PR_ALL, TYPE1, OTUS, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

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

def uniprot_label(h, tag):
    # >sp|Q9F7P4|PRRG_PRB01 Green-light absorbing proteorhodopsin OS=... OX=...
    m = re.match(r'(sp|tr)\|([^|]+)\|(\S+)\s+(.*)', h)
    if m:
        acc, rest = m.group(2), m.group(4)
        osm = re.search(r'OS=([^O]+?)(\s+OX=|$)', rest)
        osname = osm.group(1).strip()[:28].replace(' ', '_') if osm else 'env'
        return '%s|%s|%s' % (acc, tag, osname)
    return '%s|%s|env' % (h.split()[0][:15], tag)

T1_CLASS = [
    (r'xanthorhodopsin', 'xR'), (r'bacteriorhodopsin', 'bR'),
    (r'halorhodopsin', 'HR'), (r'sensory rhodopsin', 'sR'),
    (r'actinorhodopsin', 'ActR'), (r'channelrhodopsin', 'ChR'),
    (r'archaerhodopsin', 'ArchR'), (r'deltarhodopsin', 'dR'),
    (r'cruxrhodopsin', 'cR'),
]

out_h = []
po, ps = read_fasta(PR_ALL)
for h in po:
    hl = h.lower()
    if 'green' in hl and 'proteorhodopsin' in hl: tag = 'GPR'
    elif 'blue' in hl and 'proteorhodopsin' in hl: tag = 'BPR'
    else: tag = 'PR'
    out_h.append((uniprot_label(h, tag), ps[h]))

to, ts = read_fasta(TYPE1)
for h in to:
    hl = h.lower()
    tag = 'Rho'
    for pat, t in T1_CLASS:
        if re.search(pat, hl): tag = t; break
    out_h.append((uniprot_label(h, tag), ts[h]))

oo, os_ = read_fasta(OTUS)
for h in oo:
    out_h.append((h, os_[h]))

# write
with open(OUT, 'w') as f:
    for label, seq in out_h:
        label = label.replace(' ', '_')
        f.write('>%s\n' % label)
        for j in range(0, len(seq), 60): f.write(seq[j:j+60] + '\n')

tags = Counter(l.split('|')[1] for l, _ in out_h)
print('[done] tree input: %d sequences' % len(out_h))
print('  label breakdown: %s' % dict(tags))
