#!/usr/bin/env python3
"""How many clusters of placed SFS (w >= 5, full germline set, clonal and subclonal regions) of each somatic class are
the source of a candidate SV: S tokens of the CALL rows of R/bp_coverage/nogl_w5/provenance.tsv matched to the
intervals of clusters.full.tsv.gz (sfs_classes.py). Each cluster counts once.
Output: cluster_to_sv.tsv. Usage: python cluster_to_sv.py"""
import gzip, os, re
from collections import Counter
O = os.environ["HG008"] + "/reads_ec_align/sfs_signal"; os.makedirs(O, exist_ok=True)
PROV = os.environ["HG008"] + "/reads_ec_align/bp_coverage/nogl_w5/provenance.tsv"
cl = {}
for l in gzip.open(f"{O}/clusters.full.tsv.gz", "rt"):
    if l.startswith("chrom"): continue
    c, s, e, w, k = l.rstrip("\n").split("\t")
    if int(w) >= 5: cl[(c, int(s), int(e))] = k
src = set()
for l in open(PROV):
    f = l.rstrip("\n").split("\t")
    if f[0] != "CALL": continue
    for m in re.finditer(r"S:([^:|]+):(\d+)-(\d+)", f[8]):
        key = (m.group(1), int(m.group(2)), int(m.group(3)) + 1)
        if key in cl: src.add(key)
tot = Counter(cl.values()); sv = Counter(cl[k] for k in src)
with open(f"{O}/cluster_to_sv.tsv", "w") as fo:
    fo.write("class\tclusters\tsource_of_candidate_sv\n")
    for k in ("SNV", "INDEL_HP", "INDEL_other", "no_smvar", "other_trf", "other"):
        fo.write(f"{k}\t{tot[k]}\t{sv[k]}\n")
