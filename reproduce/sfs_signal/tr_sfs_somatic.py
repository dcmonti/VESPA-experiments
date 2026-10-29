#!/usr/bin/env python3
"""For each breakpoint without SFS (lost_points.tsv), the SFS of the tumor reads that carry the somatic allele.

Somatic allele: a length (delta of tr_alleles.tsv, error-corrected reads) carried by at least MIN_READS tumor reads and
by no normal read within +- TOL bp. For each tumor read with that allele, its placed SFS (sfs_clean.bed.gz of the index:
the SFS the pipeline uses, placement through the tumor EC alignments) that overlap the locus of tr_alleles.py, and the
distance of the nearest of them from the breakpoint (0 = overlapping it).
Indexes: the full germline set (R/bp_coverage) and reference + HG008-N v6.3 (R/germline_cascade/rows/R2w_refv63_nofilters).
Output: tr_sfs_somatic.tsv (one row per breakpoint and index). Usage: python tr_sfs_somatic.py
"""
import csv, os, subprocess
from collections import Counter, defaultdict

MIN_READS, TOL = 3, 3
R = os.environ["HG008"] + "/reads_ec_align"
O = os.environ["HG008"] + "/reads_ec_align/sfs_signal"; os.makedirs(O, exist_ok=True)
IDX = [("full", f"{R}/bp_coverage/sfs_clean.bed.gz"),
       ("v63", f"{R}/germline_cascade/rows/R2w_refv63_nofilters/sfs_clean.bed.gz")]

pts = [r for r in csv.DictReader(open(f"{O}/lost_points.tsv"), delimiter="\t")]
ids = {r["id"] for r in pts}
rows = [r for r in csv.DictReader(open(f"{O}/tr_alleles.tsv"), delimiter="\t") if r["reads"] == "ec" and r["id"] in ids]
som, locus = {}, {}
for id_ in ids:
    t = Counter(int(r["delta"]) for r in rows if r["id"] == id_ and r["sample"] == "tumor")
    n = [int(r["delta"]) for r in rows if r["id"] == id_ and r["sample"] == "normal"]
    som[id_] = sorted(d for d, c in t.items() if c >= MIN_READS and all(abs(d - x) > TOL for x in n))
    r0 = next(r for r in rows if r["id"] == id_); locus[id_] = (r0["chrom"], int(r0["locus_start"]), int(r0["locus_end"]))
reads = defaultdict(set)  # read -> ids
for r in rows:
    if r["sample"] == "tumor" and int(r["delta"]) in som[r["id"]]: reads[r["read"]].add(r["id"])
nsom = {id_: len({r["read"] for r in rows if r["id"] == id_ and r["sample"] == "tumor" and int(r["delta"]) in som[id_]})
        for id_ in ids}
with open(f"{O}/tr_sfs_somatic.tsv", "w") as fo:
    fo.write("index\tid\trole\tchrom\tpos\tsomatic_alleles\treads_with_allele\treads_with_sfs\tnearest_sfs_bp\n")
    for name, path in IDX:
        plc = defaultdict(list)
        zc = subprocess.Popen(["zcat", path], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for l in zc.stdout:
            f = l.split("\t", 5)
            if f[3] in reads:
                for id_ in reads[f[3]]:
                    ch, s, e = locus[id_]
                    if f[0] == ch and int(f[1]) < e and int(f[2]) > s: plc[id_].append((f[3], int(f[1]), int(f[2])))
        zc.wait()
        for r in pts:
            id_, p = r["id"], int(r["pos"]) - 1
            hits = plc.get(id_, [])
            near = min((max(s - p, p - e + 1, 0) for _, s, e in hits), default="-")
            fo.write(f"{name}\t{id_}\t{r['role']}\t{r['chrom']}\t{r['pos']}\t{','.join(map(str, som[id_])) or '-'}\t"
                     f"{nsom[id_]}\t{len({rd for rd, _, _ in hits})}\t{near}\n")
print("done")
