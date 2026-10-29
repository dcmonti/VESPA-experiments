#!/usr/bin/env python3
"""Homopolymer indels of smvar_cov.tsv.gz (sfs_classes.py): length of the homopolymer run in each allele. With runs
capped at CAP by homopolymer compression, the two alleles are the same string when both runs are at least CAP long.
Output: hp_runs.tsv (per variant) and a summary on stdout, with the SFS coverage (full germline set) of each group.
Usage: python hp_runs.py [CAP=2]"""
import gzip, os, sys
import pysam
CAP = int(sys.argv[1]) if len(sys.argv) > 1 else 2
O = os.environ["HG008"] + "/reads_ec_align/sfs_signal"; os.makedirs(O, exist_ok=True)
fa = pysam.FastaFile(os.environ["HG008"] + "/reads_ec_align/ref/GRCh38_GIABv3.chr1-22X.fa")
n = {True: 0, False: 0}; cov = {True: 0, False: 0}
with gzip.open(f"{O}/smvar_cov.tsv.gz", "rt") as fh, open(f"{O}/hp_runs.tsv", "w") as fo:
    h = next(fh).rstrip("\n").split("\t"); fo.write("chrom\tpos\tref\talt\tref_run\talt_run\tsame_after_hpc\tfull_sfs\n")
    for l in fh:
        r = dict(zip(h, l.rstrip("\n").split("\t")))
        if r["type"] != "INDEL_HP": continue
        p = int(r["pos"]) - 1; ref, alt = r["ref"], r["alt"]
        d = alt[1:] if len(alt) > len(ref) else ref[1:]; b = d[0].upper()
        # run of b in the reference around the indel (after the anchor base, extended both ways)
        s = e = p + 1
        while s > 0 and fa.fetch(r["chrom"], s - 1, s).upper() == b: s -= 1
        while fa.fetch(r["chrom"], e, e + 1).upper() == b: e += 1
        rr = e - s; ar = rr + len(d) if len(alt) > len(ref) else rr - len(d)
        same = rr >= CAP and ar >= CAP
        n[same] += 1; cov[same] += r["full_sfs"] == "1"
        fo.write(f"{r['chrom']}\t{r['pos']}\t{ref}\t{alt}\t{rr}\t{ar}\t{int(same)}\t{r['full_sfs']}\n")
print(f"CAP={CAP}: same after compression {n[True]} (with SFS {cov[True]}); different {n[False]} (with SFS {cov[False]})")
