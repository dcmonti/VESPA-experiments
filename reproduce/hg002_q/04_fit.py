#!/usr/bin/env python3
"""HG002, chr1: the per-read detection rate q of LOWPOWER, from visibility.chr1.tsv (03_visibility.py).

For each heterozygous SV, D = spanning reads tagged on the carrying homolog and K = those that show it (strict and
loose criteria). Only tagged reads enter, since the homolog of an untagged read is unknown.
  q        pooled K / D, with a bootstrap 95% interval over the SVs; by type, length class and tandem repeat
  check 1  independence of the reads (assumption (ii)): the fraction of SVs with K = 0 against the binomial
           prediction mean (1 - q)^D, by class of D; and a beta-binomial fit of (K, D), whose intra-locus
           correlation rho = 1 / (a + b + 1) is 0 when the reads of a locus fail independently (likelihood-ratio test
           against the binomial)
  check 2  the other homolog (assumption (i)): reads tagged on the homolog that does not carry the SV and show it
  check 3  homozygous SVs: K / D over all spanning reads, tagged or not
Output: q_summary.tsv, q_by_depth.tsv. Usage: python 04_fit.py OUT_DIR
"""
import csv, sys
import numpy as np
from scipy import optimize, stats, special

OUT = sys.argv[1]
CHROM = "chr1"
rng = np.random.default_rng(20261016)
rows = list(csv.DictReader(open(f"{OUT}/visibility.{CHROM}.tsv"), delimiter="\t"))
for r in rows:
    for k in r:
        if k.startswith(("span_", "strict_", "loose_")) or k in ("len", "in_trf"):
            r[k] = int(r[k])
het = [r for r in rows if r["carrier"] in ("1", "2")]
hom = [r for r in rows if r["carrier"] == "both"]


def carrier_counts(r, crit):
    h, o = f"hp{r['carrier']}", f"hp{3 - int(r['carrier'])}"
    return r[f"span_{h}"], r[f"{crit}_{h}"], r[f"span_{o}"], r[f"{crit}_{o}"]


def pooled(D, K, boot=1000):
    D, K = np.asarray(D), np.asarray(K)
    if D.sum() == 0:
        return np.nan, np.nan, np.nan
    q = K.sum() / D.sum()
    idx = rng.integers(0, len(D), (boot, len(D)))
    qs = K[idx].sum(1) / np.maximum(D[idx].sum(1), 1)
    return q, *np.quantile(qs, [0.025, 0.975])


def betabinom_fit(D, K):
    """MLE of Beta(a, b) for the locus rate; returns a, b, rho, LR statistic and p against the binomial."""
    D, K = np.asarray(D), np.asarray(K)
    def nll(t):
        a, b = np.exp(t)
        return -stats.betabinom.logpmf(K, D, a, b).sum()
    q0 = K.sum() / D.sum()
    res = optimize.minimize(nll, np.log([q0 * 10, (1 - q0) * 10]), method="Nelder-Mead")
    a, b = np.exp(res.x)
    ll_bin = stats.binom.logpmf(K, D, q0).sum()
    lr = 2 * (-res.fun - ll_bin)
    return a, b, 1 / (a + b + 1), lr, 0.5 * stats.chi2.sf(lr, 1)   # boundary: half chi2(1)


summ = open(f"{OUT}/q_summary.tsv", "w")
summ.write("criterion\tsubset\tn_sv\treads\tq\tq_lo\tq_hi\tother_homolog_rate\n")
bydep = open(f"{OUT}/q_by_depth.tsv", "w")
bydep.write("criterion\tD\tn_sv\tobs_frac_K0\tpred_binomial\tpred_betabinomial\n")
for crit in ("strict", "loose"):
    subsets = {"all": het,
               "DEL": [r for r in het if r["type"] == "DEL"], "INS": [r for r in het if r["type"] == "INS"],
               "50-300 bp": [r for r in het if r["len"] < 300], "300-1000 bp": [r for r in het if 300 <= r["len"] < 1000],
               "1-10 kb": [r for r in het if r["len"] >= 1000],
               "tandem repeat": [r for r in het if r["in_trf"]], "no tandem repeat": [r for r in het if not r["in_trf"]]}
    for name, sub in subsets.items():
        c = [carrier_counts(r, crit) for r in sub]
        c = [x for x in c if x[0] > 0]
        if not c:
            continue
        D, K, Do, Ko = map(np.array, zip(*c))
        q, lo, hi = pooled(D, K)
        other = Ko.sum() / Do.sum() if Do.sum() else np.nan
        summ.write(f"{crit}\t{name}\t{len(c)}\t{D.sum()}\t{q:.4f}\t{lo:.4f}\t{hi:.4f}\t{other:.4f}\n")
        if name != "all":
            continue
        a, b, rho, lr, p = betabinom_fit(D, K)
        summ.write(f"{crit}\tbeta-binomial a={a:.3f} b={b:.3f} rho={rho:.4f} LR={lr:.2f} p={p:.3g}\t\t\t\t\t\t\n")
        for d in sorted(set(D)):
            m = D == d
            if m.sum() < 5:
                continue
            obs = (K[m] == 0).mean()
            bb0 = special.beta(a, b + d) / special.beta(a, b)
            bydep.write(f"{crit}\t{d}\t{m.sum()}\t{obs:.4f}\t{(1 - q) ** d:.4g}\t{bb0:.4g}\n")
    # homozygous: all spanning reads
    D = np.array([r["span_hp1"] + r["span_hp2"] + r["span_hp0"] for r in hom])
    K = np.array([r[f"{crit}_hp1"] + r[f"{crit}_hp2"] + r[f"{crit}_hp0"] for r in hom])
    m = D > 0
    if m.any():
        q, lo, hi = pooled(D[m], K[m])
        summ.write(f"{crit}\thomozygous, all reads\t{m.sum()}\t{D[m].sum()}\t{q:.4f}\t{lo:.4f}\t{hi:.4f}\t\n")
summ.close(); bydep.close()
print(open(f"{OUT}/q_summary.tsv").read())
