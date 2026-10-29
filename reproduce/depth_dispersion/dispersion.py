#!/usr/bin/env python3
"""Over-dispersion of the read depth of a HiFi alignment, under the model of the error-correction threshold.

The depth X at a locus is modelled as a negative binomial with mean mu and variance mu + alpha * mu^2.
This script measures alpha on a set of random autosomal positions:
  1. positions: drawn uniformly on chr1-22 (fixed seed), rejecting N bases of the reference and the
     centromere/telomere intervals excluded by VESPA call; the same positions are used for every sample
  2. depth: number of alignments overlapping each position, excluding unmapped, secondary, QC-fail and
     duplicate records (the default of samtools depth); a read deleting the position still covers it
  3. alpha, on all positions and on the positions outside tandem repeats:
       moments   (var - mean) / mean^2
       ML        maximum likelihood of NB(mean, alpha), mean fixed at the sample mean (its ML estimate)
     each with a bootstrap 95% interval over positions

Usage:
  dispersion.py positions --ref REF.fa.gz --exclude BED --n N --seed S --out positions.tsv.gz
  dispersion.py depth --bam BAM --positions positions.tsv.gz --threads T --out depth.tsv.gz
  dispersion.py fit --depth depth.tsv.gz --trf BED --sample NAME --boot B --out summary.tsv
"""
import argparse
import bisect
import gzip
import sys
from multiprocessing import Pool

import numpy as np
import pysam
from scipy import optimize, special

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def read_bed(path, chroms):
    iv = {c: [] for c in chroms}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if f[0] in iv:
                iv[f[0]].append((int(f[1]), int(f[2])))
    out = {}
    for c, lst in iv.items():
        lst.sort()
        merged = []
        for s, e in lst:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[c] = ([s for s, _ in merged], [e for _, e in merged])
    return out


def inside(iv, chrom, pos):
    starts, ends = iv[chrom]
    i = bisect.bisect_right(starts, pos) - 1
    return i >= 0 and pos < ends[i]


def cmd_positions(a):
    fa = pysam.FastaFile(a.ref)
    lengths = np.array([fa.get_reference_length(c) for c in AUTOSOMES], dtype=np.int64)
    excl = read_bed(a.exclude, AUTOSOMES)
    rng = np.random.default_rng(a.seed)
    kept = set()
    while len(kept) < a.n:
        k = 2 * (a.n - len(kept))
        idx = rng.choice(len(AUTOSOMES), size=k, p=lengths / lengths.sum())
        pos = (rng.random(k) * lengths[idx]).astype(np.int64)
        for i, p in zip(idx, pos):
            c = AUTOSOMES[i]
            if (c, p) in kept or inside(excl, c, p):
                continue
            if fa.fetch(c, int(p), int(p) + 1).upper() == "N":
                continue
            kept.add((c, int(p)))
            if len(kept) == a.n:
                break
    rows = sorted(kept, key=lambda x: (AUTOSOMES.index(x[0]), x[1]))
    with gzip.open(a.out, "wt") as out:
        out.write("chrom\tpos0\n")
        for c, p in rows:
            out.write(f"{c}\t{p}\n")


def depth_chunk(args):
    bam_path, chrom, positions = args
    bam = pysam.AlignmentFile(bam_path)
    # read_callback "all" skips unmapped, secondary, QC-fail and duplicate records
    return [(chrom, p, bam.count(chrom, p, p + 1, read_callback="all")) for p in positions]


def cmd_depth(a):
    by_chrom = {}
    with gzip.open(a.positions, "rt") as fh:
        next(fh)
        for line in fh:
            c, p = line.split()
            by_chrom.setdefault(c, []).append(int(p))
    jobs = []
    for c, ps in by_chrom.items():
        for i in range(0, len(ps), 1000):
            jobs.append((a.bam, c, ps[i:i + 1000]))
    with Pool(a.threads) as pool:
        res = pool.map(depth_chunk, jobs)
    with gzip.open(a.out, "wt") as out:
        out.write("chrom\tpos0\tdepth\n")
        for chunk in res:
            for c, p, d in chunk:
                out.write(f"{c}\t{p}\t{d}\n")


def alpha_moments(x):
    m, v = x.mean(), x.var(ddof=1)
    return (v - m) / m ** 2


def nb_loglik(x, mu, alpha):
    r = 1.0 / alpha
    return np.sum(special.gammaln(x + r) - special.gammaln(r) - special.gammaln(x + 1)
                  + r * np.log(r / (r + mu)) + x * np.log(mu / (r + mu)))


def alpha_ml(x):
    mu = x.mean()
    res = optimize.minimize_scalar(lambda la: -nb_loglik(x, mu, np.exp(la)),
                                   bounds=(np.log(1e-6), np.log(10.0)), method="bounded")
    return float(np.exp(res.x))


def cmd_fit(a):
    trf = read_bed(a.trf, AUTOSOMES)
    chrom, depth, in_trf = [], [], []
    with gzip.open(a.depth, "rt") as fh:
        next(fh)
        for line in fh:
            c, p, d = line.split()
            depth.append(int(d))
            in_trf.append(inside(trf, c, int(p)))
    depth = np.array(depth, dtype=np.float64)
    in_trf = np.array(in_trf)
    rng = np.random.default_rng(1)
    cols = ["sample", "positions", "n", "mean", "variance", "median",
            "alpha_moments", "alpha_moments_lo", "alpha_moments_hi",
            "alpha_ml", "alpha_ml_lo", "alpha_ml_hi"]
    with open(a.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for label, x in [("all", depth), ("outside_trf", depth[~in_trf])]:
            bm, bl = [], []
            for _ in range(a.boot):
                s = x[rng.integers(0, len(x), len(x))]
                bm.append(alpha_moments(s))
                bl.append(alpha_ml(s))
            row = [a.sample, label, len(x), f"{x.mean():.3f}", f"{x.var(ddof=1):.3f}", f"{np.median(x):.0f}",
                   f"{alpha_moments(x):.5f}", *(f"{q:.5f}" for q in np.percentile(bm, [2.5, 97.5])),
                   f"{alpha_ml(x):.5f}", *(f"{q:.5f}" for q in np.percentile(bl, [2.5, 97.5]))]
            out.write("\t".join(map(str, row)) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("positions")
    p.add_argument("--ref", required=True)
    p.add_argument("--exclude", required=True)
    p.add_argument("--n", type=int, default=100000)
    p.add_argument("--seed", type=int, default=20261012)
    p.add_argument("--out", required=True)
    p = sub.add_parser("depth")
    p.add_argument("--bam", required=True)
    p.add_argument("--positions", required=True)
    p.add_argument("--threads", type=int, default=22)
    p.add_argument("--out", required=True)
    p = sub.add_parser("fit")
    p.add_argument("--depth", required=True)
    p.add_argument("--trf", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--boot", type=int, default=200)
    p.add_argument("--out", required=True)
    a = ap.parse_args()
    {"positions": cmd_positions, "depth": cmd_depth, "fit": cmd_fit}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
