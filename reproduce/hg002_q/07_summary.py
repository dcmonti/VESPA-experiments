#!/usr/bin/env python3
"""HG002: summary of the detection rate q over several chromosomes (the per-chromosome files of 03 and 06).

For the heterozygous SVs, against the representation in the reads (visibility_reads.<chrom>.tsv, strict and loose
criterion) and against the truth record (visibility.<chrom>.tsv, strict): q = K / D on the carrying homolog, the rate
of the other homolog, and the miss probability at d reads drawn without replacement at each locus, with the upper end
of a bootstrap 95% interval over the loci; for the read representation also the beta-binomial intra-locus correlation
rho. Subsets: all, isolated (no other truth variant >= 10 bp within 500 bp), with nearby variants, outside / in
tandem repeats, and the GC-rich chromosomes 16, 17, 19, 22. Homozygous SVs: q over all spanning reads.

Output: <dir>/q_summary.tsv, printed. Usage: python 07_summary.py RESULT_DIR chr1 chr2 ...
"""
import csv, sys
import numpy as np
from scipy import optimize, stats
from scipy.special import comb

D_DIR, CHROMS = sys.argv[1], sys.argv[2:]
GC_RICH = {"chr16", "chr17", "chr19", "chr22"}
DS = (1, 2, 3, 4, 6)
rng = np.random.default_rng(20261016)


def read(name):
    rows = []
    for c in CHROMS:
        for r in csv.DictReader(open(f"{D_DIR}/{name}.{c}.tsv"), delimiter="\t"):
            r["chrom"] = c
            rows.append(r)
    return rows


def curve(D, K):
    out = []
    for d in DS:
        v = np.array([comb(a - b, d) / comb(a, d) for a, b in zip(D, K) if a >= d])
        bs = rng.choice(v, (500, len(v))).mean(1)
        out += [v.mean(), np.quantile(bs, 0.975)]
    return out


def rho(D, K):
    q0 = K.sum() / D.sum()
    nll = lambda t: -stats.betabinom.logpmf(K, D, *np.exp(t)).sum()
    a, b = np.exp(optimize.minimize(nll, np.log([q0 * 10, (1 - q0) * 10]), method="Nelder-Mead").x)
    return 1 / (a + b + 1)


reads = read("visibility_reads")
truth = [r for r in read("visibility") if r["carrier"] in ("1", "2")]
key = lambda r: (r["chrom"], r["pos"], r["type"], r["len"])
iso = {key(r): r["isolated"] == "1" for r in reads}
for r in truth:
    r["isolated"] = "1" if iso.get(key(r)) else "0"
subsets = {"all": lambda r: True, "isolated": lambda r: r["isolated"] == "1",
           "with nearby variants": lambda r: r["isolated"] == "0",
           "outside tandem repeats": lambda r: r["in_trf"] == "0", "in tandem repeats": lambda r: r["in_trf"] == "1",
           "GC-rich chromosomes": lambda r: r["chrom"] in GC_RICH}
cols = ["representation", "criterion", "subset", "n", "reads", "q", "other_homolog", "rho"] + \
       [f"miss_d{d}{s}" for d in DS for s in ("", "_hi")]
out = open(f"{D_DIR}/q_summary.tsv", "w")
out.write("\t".join(cols) + "\n")
lines = []
for rep, rows, crits in (("reads", reads, (("strict", "K_strict", "K_other_strict"), ("loose", "K_loose", "K_other_loose"))),
                         ("truth", truth, (("strict", None, None),))):
    for crit, kc, koc in crits:
        for name, f in subsets.items():
            sub = [r for r in rows if f(r)]
            if rep == "reads":
                x = np.array([[int(r["D"]), int(r[kc]), int(r["D_other"]), int(r[koc])] for r in sub if int(r["D"]) > 0])
            else:
                x = []
                for r in sub:
                    h, o = f"hp{r['carrier']}", f"hp{3 - int(r['carrier'])}"
                    if int(r[f"span_{h}"]) > 0:
                        x.append([int(r[f"span_{h}"]), int(r[f"strict_{h}"]), int(r[f"span_{o}"]), int(r[f"strict_{o}"])])
                x = np.array(x)
            if len(x) == 0:
                continue
            D, K, Do, Ko = x.T
            row = [rep, crit, name, len(x), int(D.sum()), K.sum() / D.sum(), Ko.sum() / max(Do.sum(), 1),
                   rho(D, K) if rep == "reads" and name == "all" else np.nan] + curve(D, K)
            out.write("\t".join(f"{v:.4f}" if isinstance(v, float) else str(v) for v in row) + "\n")
            lines.append(row)
hom = [r for r in read("visibility") if r["carrier"] == "both"]
D = np.array([sum(int(r[f"span_hp{h}"]) for h in (1, 2, 0)) for r in hom])
K = np.array([sum(int(r[f"strict_hp{h}"]) for h in (1, 2, 0)) for r in hom])
m = D > 0
out.write(f"truth\tstrict\thomozygous, all reads\t{m.sum()}\t{D[m].sum()}\t{K[m].sum() / D[m].sum():.4f}\n")
out.close()
print(open(f"{D_DIR}/q_summary.tsv").read())
