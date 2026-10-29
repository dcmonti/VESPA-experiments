#!/usr/bin/env python3
"""Fraction of the regions where the classes 'other, tandem repeat' and 'other' of sfs_classes.py can occur (V0.5 clonal
and subclonal regions inside the smvar v0.3 regions) that lies within 50 bp of a tandem repeat of the Severus TRF/VNTR
annotation. Output: trf_fraction.txt. Usage: python trf_fraction.py"""
import os
from collections import defaultdict
A = os.environ["HG008"]
O = os.environ["HG008"] + "/reads_ec_align/sfs_signal"; os.makedirs(O, exist_ok=True)
def bed(p, pad=0):
    d = defaultdict(list)
    for l in open(p):
        f = l.split("\t"); d[f[0]].append((max(0, int(f[1]) - pad), int(f[2]) + pad))
    out = {}
    for c, L in d.items():
        L.sort(); m = []
        for s, e in L:
            if m and s <= m[-1][1]: m[-1][1] = max(m[-1][1], e)
            else: m.append([s, e])
        out[c] = m
    return out
def inter(a, b):
    out = {}
    for c in set(a) & set(b):
        x, y = a[c], b[c]; i = j = 0; m = []
        while i < len(x) and j < len(y):
            s, e = max(x[i][0], y[j][0]), min(x[i][1], y[j][1])
            if s < e: m.append([s, e])
            if x[i][1] < y[j][1]: i += 1
            else: j += 1
        out[c] = m
    return out
size = lambda d: sum(e - s for v in d.values() for s, e in v)
reg = bed(f"{A}/ground_truths/v0.5/GRCh38_HG008-T-V0.5_somatic-stvar-clonal_and_subclonal.draftbenchmark.bed")
sm = bed(f"{A}/ground_truths/smvar_v0.3/HG008-T_somatic_smvar_benchmark_v0.3_all.bed")
trf = bed(f"{A}/refs/sev_extras/human_GRCh38_no_alt_analysis_set.trf.bed", 50)
both = inter(reg, sm); t = size(inter(both, trf))
with open(f"{O}/trf_fraction.txt", "w") as fo:
    fo.write(f"regions {size(reg)} bp; inside smvar regions {size(both)} bp; within 50 bp of a tandem repeat {t} bp "
             f"({t / size(both):.4f})\n")
