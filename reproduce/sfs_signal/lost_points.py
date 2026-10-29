#!/usr/bin/env python3
"""The breakpoints of the V0.5 truth (clonal and subclonal regions) that no SFS of the full germline set covers, and
whether the tumor and normal EC reads carry an allele of the event there.

A read carries the allele when, among the primary alignments with MAPQ >= 20 over the breakpoint, it has an indel at
least half as long as the SV within 100 bp + the SV length of the breakpoint, or a soft clip within 50 bp of it.
The test accepts any indel of similar length, so in a tandem repeat it also counts germline length alleles.
For the SFS of the reference-only index that covered the breakpoint, the fraction of their reads that carry the allele
tells whether that coverage came from the event or from germline differences nearby.

Input: bp_points.tsv, sfs_hits.tsv.gz (bp_points.py). Output: lost_points.tsv
Usage: python lost_points.py
"""
import csv, gzip, os
from collections import defaultdict
import pysam

O = os.environ["HG008"] + "/reads_ec_align/sfs_signal"; os.makedirs(O, exist_ok=True)
R = os.environ["HG008"] + "/reads_ec_align/align"


def carriers(bam, ch, p, L, names=None):
    tot = car = sn = sc = 0
    for a in bam.fetch(ch, p - 1, p):
        if a.is_secondary or a.is_supplementary or a.mapping_quality < 20: continue
        tot += 1; rp = a.reference_start; c = False; ct = a.cigartuples
        if ct[0][0] == 4 and abs(a.reference_start - p) <= 50: c = True
        if ct[-1][0] == 4 and abs(a.reference_end - p) <= 50: c = True
        for op, n in ct:
            if op in (0, 7, 8): rp += n
            elif op == 1 and n >= L / 2 and abs(rp - p) <= 100 + L: c = True
            elif op == 2:
                if n >= L / 2 and rp <= p + 100 + L and rp + n >= p - 100 - L: c = True
                rp += n
        car += c
        if names is not None and a.query_name in names: sn += 1; sc += c
    return tot, car, sn, sc


if __name__ == "__main__":
    rows = [r for r in csv.DictReader(open(f"{O}/bp_points.tsv"), delimiter="\t") if r["in_cs"] == "1"]
    lost = [r for r in rows if int(r["full_sfs"]) == 0]
    keys = {f"{r['id']}:{r['role']}" for r in lost}; reads = defaultdict(set)
    for l in gzip.open(f"{O}/sfs_hits.tsv.gz", "rt"):
        f = l.rstrip("\n").split("\t")
        if f[0] == "ref" and f[1] in keys: reads[f[1]].add(f[6])
    T = pysam.AlignmentFile(f"{R}/tumor.ec.hp.bam"); N = pysam.AlignmentFile(f"{R}/normal.ec.hp.bam")
    with open(f"{O}/lost_points.tsv", "w") as fo:
        fo.write("id\trole\tchrom\tpos\tsvtype\teventtype\tsvlen\tvaf\tsubclonal\tin_clonal\tin_trf\tref_sfs\t"
                 "tumor_depth\ttumor_carriers\tref_sfs_reads_here\tref_sfs_reads_carrying\tnormal_depth\tnormal_carriers\n")
        for r in lost:
            ch, p, L = r["chrom"], int(r["pos"]), int(r["svlen"] or 50); k = f"{r['id']}:{r['role']}"
            tt, tc, sn, sc = carriers(T, ch, p, L, reads[k]); nt, nc, _, _ = carriers(N, ch, p, L)
            fo.write("\t".join([r["id"], r["role"], ch, r["pos"], r["svtype"], r["eventtype"], r["svlen"], r["vaf"],
                                r["subclonal"], r["in_clonal"], r["in_trf"], r["ref_sfs"]] + list(map(str, (tt, tc, sn, sc, nt, nc)))) + "\n")
    print("done")
