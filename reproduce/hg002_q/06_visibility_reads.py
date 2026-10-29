#!/usr/bin/env python3
"""HG002, chr1: visibility of germline SVs against their REPRESENTATION IN THE READS rather than in the truth.

03_visibility.py compares each read with the truth record. Where dipcall writes one event as several records (an INS
and a DEL 2 bp apart for a complex allele) or places it differently in a repeat, no read matches the truth although
every read shows the event. In the method, the call of a germline event is built from reads aligned like the normal
ones, so the representation that the normal reads are tested against is that of the reads. Here, for each SV of the
truth (its locus and its carrying homolog only):
  representation, leave-one-out: for carrier read i, the other carrier reads that cover the breakpoint each give their
      longest indel (>= 10 bp) starting in [s - 500, s + len + 500]; the representation is the most frequent type among
      them and the median length of that type. Fewer than 2 such reads: the locus shows no event in the reads, and i
      counts as not showing it.
  read i shows the event when it has an op of that type in the window with a length that passes the strict criterion
      (min/max >= 0.90, or |diff| < 50 with length >= 10) or the loose one (|diff| <= max(50, 0.25 len)).
  other homolog: each read of the non-carrying homolog against the representation of all the carrier reads.
Reads as in 03_visibility.py (primary, MAPQ >= 60, covering the breakpoint, truth-phased HP).

Output: visibility_reads.chr1.tsv (per SV: carrier D, K strict/loose, other D, K strict/loose) and a printed summary
with the subsampled miss probability at d reads. Usage: python 06_visibility_reads.py OUT_DIR [CHROM [TAG]]
"""
import bisect, csv, sys
from collections import Counter
import numpy as np, pysam
from scipy.special import comb

OUT = sys.argv[1]
CHROM = sys.argv[2] if len(sys.argv) > 2 else "chr1"
TAG = sys.argv[3] if len(sys.argv) > 3 else CHROM
MIN_MAPQ, WIN, MIN_OP = 60, 500, 10
BAM = f"{OUT}/align/hg002.{TAG}.raw.hp.bam"
OUTD = OUT if TAG == CHROM else f"{OUT}/{TAG}"   # per-chromosome results of a multi-chromosome run
DIP = f"{OUT}/truth/GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.vcf.gz"


def indels(aln, lo, hi):
    """(type, length) of the indel ops >= MIN_OP bp starting in [lo, hi]; type 1 = INS, 2 = DEL."""
    out, rpos = [], aln.reference_start
    for op, n in aln.cigartuples:
        if op in (1, 2) and n >= MIN_OP and lo <= rpos <= hi:
            out.append((op, n))
        if op in (0, 2, 3, 7, 8):
            rpos += n
        if rpos > hi:
            break
    return out


def longest(ops):
    return max(ops, key=lambda x: x[1]) if ops else None


def representation(longs):
    longs = [x for x in longs if x]
    if len(longs) < 2:
        return None
    t = Counter(x[0] for x in longs).most_common(1)[0][0]
    return t, float(np.median([x[1] for x in longs if x[0] == t]))


def shows(ops, rep):
    if rep is None:
        return False, False
    t, L = rep
    st = lo = False
    for op, n in ops:
        if op != t:
            continue
        st = st or (n >= 50 and min(n, L) / max(n, L) >= 0.90) or (n >= 10 and abs(L - n) < 50)
        lo = lo or (n >= 50 and abs(n - L) <= max(50, 0.25 * L))
    return st, lo


def main():
    rows = list(csv.DictReader(open(f"{OUTD}/visibility.{CHROM}.tsv"), delimiter="\t"))
    ind = sorted(rec.pos for rec in pysam.VariantFile(DIP).fetch(CHROM)
                 if any(abs(len(rec.ref) - len(a)) >= 10 for a in (rec.alts or ())))
    isolated = lambda p: (bisect.bisect_right(ind, p + 500) - bisect.bisect_left(ind, p - 500)) <= 1
    bam = pysam.AlignmentFile(BAM)
    out = open(f"{OUTD}/visibility_reads.{CHROM}.tsv", "w")
    out.write("pos\ttype\tlen\tcarrier\tin_trf\tisolated\tD\tK_strict\tK_loose\tD_other\tK_other_strict\tK_other_loose\n")
    res = []
    for r in rows:
        if r["carrier"] not in ("1", "2"):
            continue
        s, L, c = int(r["pos"]), int(r["len"]), int(r["carrier"])
        lo, hi = s - WIN, s + L + WIN
        car, oth = [], []
        for aln in bam.fetch(CHROM, s, s + 1):
            if aln.is_secondary or aln.is_supplementary or aln.is_unmapped or aln.mapping_quality < MIN_MAPQ:
                continue
            if aln.reference_start > s or aln.reference_end - 1 < s or not aln.has_tag("HP"):
                continue
            ops = indels(aln, lo, hi)
            (car if aln.get_tag("HP") == c else oth).append(ops)
        longs = [longest(o) for o in car]
        k = np.array([shows(car[i], representation(longs[:i] + longs[i + 1:])) for i in range(len(car))]).reshape(-1, 2)
        rep_all = representation(longs)
        ko = np.array([shows(o, rep_all) for o in oth]).reshape(-1, 2)
        row = (s, r["type"], L, c, int(r["in_trf"]), int(isolated(s)), len(car), int(k[:, 0].sum()), int(k[:, 1].sum()),
               len(oth), int(ko[:, 0].sum()), int(ko[:, 1].sum()))
        res.append(row)
        out.write("\t".join(map(str, row)) + "\n")
    out.close()

    res = np.array([x[4:] for x in res], dtype=float)   # in_trf isolated D Ks Kl Do Kos Kol
    rng = np.random.default_rng(3)
    for name, m in (("all", np.ones(len(res), bool)), ("isolated", res[:, 1] == 1), ("non-isolated", res[:, 1] == 0),
                    ("no TR", res[:, 0] == 0), ("TR", res[:, 0] == 1)):
        x = res[m & (res[:, 2] > 0)]
        for crit, j in (("strict", 3), ("loose", 4)):
            D, K, Do, Ko = x[:, 2], x[:, j], x[:, 5], x[:, j + 3]
            curve = []
            for d in (1, 2, 3, 4, 6):
                v = np.array([comb(a - b, d) / comb(a, d) for a, b in zip(D, K) if a >= d])
                bs = [rng.choice(v, len(v)).mean() for _ in range(500)]
                curve.append(f"d={d}: {v.mean():.3f} (<={np.quantile(bs, .975):.3f})")
            print(f"{name:12s} {crit:6s} n={len(x):3d} q={K.sum() / D.sum():.3f} other={Ko.sum() / max(Do.sum(), 1):.3f}  " + "  ".join(curve))


if __name__ == "__main__":
    main()
