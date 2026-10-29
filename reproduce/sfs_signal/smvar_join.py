#!/usr/bin/env python3
"""SFS/cluster coverage (full germline set, smvar_cov.tsv.gz) of the small variants that are in our tumor reads
(>= 3 alternative reads, smvar_pileup.tsv.gz), and germline-like variants (normal depth >= 10, normal VAF >= 0.35).
Output: smvar_join_summary.tsv (and, as comment lines, the normal support of the indels the normal reads uncover).
Usage: python smvar_join.py"""
import gzip, os
O = os.environ["HG008"] + "/reads_ec_align/sfs_signal"; os.makedirs(O, exist_ok=True)
cov = {}
for l in gzip.open(f"{O}/smvar_cov.tsv.gz", "rt"):
    f = l.rstrip("\n").split("\t")
    if f[0] == "chrom": h = f; continue
    cov[tuple(f[:4])] = dict(zip(h, f))
P = [l.rstrip("\n").split("\t") for l in gzip.open(f"{O}/smvar_pileup.tsv.gz", "rt") if not l.startswith("chrom")]
with open(f"{O}/smvar_join_summary.tsv", "w") as fo:
    fo.write("type\tn\tin_reads\tin_reads_frac\tsfs_in_reads\tcluster_in_reads\tnormal_dp10\tnormal_vaf35\tnormal_vaf35_frac\n")
    for t in ("SNV", "INDEL_HP", "INDEL_other"):
        sel = [f for f in P if f[4] == t and int(f[5]) >= 0]
        inr = [f for f in sel if int(f[6]) >= 3]; c = [cov[tuple(f[:4])] for f in inr]
        s = sum(x["full_sfs"] == "1" for x in c); k = sum(x["full_cluster"] == "1" for x in c)
        nn = [f for f in sel if int(f[7]) >= 10]; g = sum(1 for f in nn if int(f[8]) / int(f[7]) >= 0.35)
        fo.write(f"{t}\t{len(sel)}\t{len(inr)}\t{len(inr)/len(sel):.4f}\t{s/len(inr):.4f}\t{k/len(inr):.4f}\t{len(nn)}\t{g}\t{g/len(nn):.4f}\n")
# Indels outside homopolymers covered by SFS with the contigs alone (utg) but not with the full germline set (the normal
# reads added), against those covered by both: fraction with at least one alternative read in the normal.
pl = {tuple(f[:4]): f for f in P}
lost = [k for k, v in cov.items() if v["type"] == "INDEL_other" and v["utg_sfs"] == "1" and v["full_sfs"] == "0"]
kept = [k for k, v in cov.items() if v["type"] == "INDEL_other" and v["utg_sfs"] == "1" and v["full_sfs"] == "1"]
with open(f"{O}/smvar_join_summary.tsv", "a") as fo:
    fo.write("# INDEL_other covered with contigs only: lost with the normal reads / kept, with >= 1 normal alt read\n")
    for lab, ks in (("lost", lost), ("kept", kept)):
        a = sum(1 for k in ks if int(pl[k][8]) >= 1)
        fo.write(f"# {lab}\t{len(ks)}\t{a}\t{a / len(ks):.4f}\n")
