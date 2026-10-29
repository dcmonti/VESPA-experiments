#!/usr/bin/env python3
"""For each clonal somatic small variant: is it in our tumor reads, and is it covered by an SFS / a cluster?

Variants: GIAB HG008-T smvar draft v0.3 tumorvariants.vcf.gz, restricted to all.bed.
Types: SNV; INDEL_HP (the indel adds or removes copies of a single base that also flanks it:
a homopolymer length change, invisible after homopolymer compression); INDEL_other.
Alt reads in our data: primary alignments, MAPQ >= 60, in the tumor EC BAM (the one given to
the call). SNV: base at POS equals ALT. Indel: the read has an indel of the same type and
length starting within +-WIN bp of POS (tolerates left/right placement in repeats).
Coverage by the pipeline: >= 1 SFS (MAPQ 60, sfs_strat/<run>.sfs.tsv.gz) or >= 1 cluster of
weight >= 5 (<run>/call/chain_w5_clusters/calls.clusters) intersecting [POS-PAD, END+PAD].
Output: per-variant table small_variant_support.tsv.gz and a summary on stdout.
USAGE: python3 small_variant_support.py [THREADS=24] [PAD=10] [WIN=5] [MIN_ALT=3]
"""
import bisect
import gzip
import os
import sys
from collections import defaultdict
from multiprocessing import Pool

import pysam

THREADS = int(sys.argv[1]) if len(sys.argv) > 1 else 24
PAD = int(sys.argv[2]) if len(sys.argv) > 2 else 10
WIN = int(sys.argv[3]) if len(sys.argv) > 3 else 5
MIN_ALT = int(sys.argv[4]) if len(sys.argv) > 4 else 3
O = os.path.dirname(os.path.abspath(__file__))
H = os.path.dirname(O)
A = os.environ["HG008"]
BAM = f"{A}/jan2024_2026-09-20/het6/align/tumor.ec.hp_tag.bam"
REF = f"{A}/jan2024_2026-09-20/refs/GIABv3.primary_noY.fa"
SM = f"{A}/ground_truths/smvar_v0.3/HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz"
BED = f"{A}/ground_truths/smvar_v0.3/HG008-T_somatic_smvar_benchmark_v0.3_all.bed"
RUNS = ["cap2", "cap2_verkko22"]


def tree(d):
    out = {}
    for c, L in d.items():
        L.sort()
        m = []
        for s, e in L:
            if m and s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        out[c] = ([a for a, _ in m], [b for _, b in m])
    return out


def hit(t, c, s, e):
    if c not in t:
        return False
    st, en = t[c]
    i = bisect.bisect_right(st, e) - 1
    return i >= 0 and en[i] >= s


def inside(t, c, s, e):
    if c not in t:
        return False
    st, en = t[c]
    i = bisect.bisect_right(st, s) - 1
    return i >= 0 and en[i] >= e


def vtype(ref, alt, fa, c, p):
    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    # anchored indel: the changed bases follow the first (anchor) base
    d = alt[1:] if len(alt) > len(ref) else ref[1:]
    if len(set(d)) == 1:
        b = d[0]
        after = fa.fetch(c, p + len(ref), p + len(ref) + 1).upper()  # base after the REF allele
        if ref[0].upper() == b or after == b:
            return "INDEL_HP"
    return "INDEL_other"


def count_alt(bam, c, p, ref, alt, t):
    """Primary MAPQ-60 reads covering POS, and those carrying the ALT."""
    dp = na = 0
    if t == "SNV":
        for col in bam.pileup(c, p, p + 1, truncate=True, min_mapping_quality=60,
                              stepper="nofilter", ignore_orphans=False, min_base_quality=0):
            for pr in col.pileups:
                a = pr.alignment
                if a.is_secondary or a.is_supplementary:
                    continue
                dp += 1
                if not pr.is_del and not pr.is_refskip and a.query_sequence[pr.query_position].upper() == alt.upper():
                    na += 1
        return dp, na
    ins = len(alt) > len(ref)
    L = abs(len(alt) - len(ref))
    for a in bam.fetch(c, p, p + 1):
        if a.is_secondary or a.is_supplementary or a.mapping_quality < 60:
            continue
        if a.reference_start > p - WIN or a.reference_end < p + len(ref) + WIN:
            continue
        dp += 1
        rp = a.reference_start
        for op, n in a.cigartuples:
            if op in (0, 7, 8):
                rp += n
            elif op == 1:
                if ins and n == L and abs(rp - (p + 1)) <= WIN:
                    na += 1
                    break
            elif op == 2:
                if not ins and n == L and abs(rp - (p + 1)) <= WIN:
                    na += 1
                    break
                rp += n
            elif op == 3:
                rp += n
            if rp > p + WIN + len(ref):
                break
    return dp, na


def work(chrom_vars):
    c, vs = chrom_vars
    out = []
    with pysam.AlignmentFile(BAM) as bam:
        for p, ref, alt, t in vs:
            dp, na = count_alt(bam, c, p, ref, alt, t)
            out.append((c, p, ref, alt, t, dp, na))
    return out


if __name__ == "__main__":
    bed = defaultdict(list)
    for line in open(BED):
        f = line.split()
        bed[f[0]].append((int(f[1]), int(f[2])))
    bed = tree(bed)
    fa = pysam.FastaFile(REF)
    byc = defaultdict(list)
    for line in gzip.open(SM, "rt"):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        p = int(f[1]) - 1
        if not inside(bed, f[0], p, p + len(f[3])):
            continue
        byc[f[0]].append((p, f[3], f[4], vtype(f[3], f[4], fa, f[0], p)))
    with Pool(THREADS) as pool:
        res = [x for r in pool.map(work, sorted(byc.items(), key=lambda kv: -len(kv[1])), chunksize=1) for x in r]

    cov = {}
    for r in RUNS:
        sf = defaultdict(list)
        with gzip.open(f"{O}/{r}.sfs.tsv.gz", "rt") as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[7] in ("trf", "nontrf", "softclip"):
                    sf[f[3]].append((int(f[4]), max(int(f[5]), int(f[4]) + 1)))
        cl = defaultdict(list)
        for line in open(f"{H}/{r}/call/chain_w5_clusters/calls.clusters"):
            f = line.split("\t", 2)
            if int(f[1]) < 5:
                continue
            c, rng = f[0].rsplit(":", 1)
            s, e = map(int, rng.split("-"))
            cl[c].append((s - 1, e))
        cov[r] = (tree(sf), tree(cl))

    with gzip.open(f"{O}/small_variant_support.tsv.gz", "wt") as fo:
        fo.write("chrom\tpos\tref\talt\ttype\tdepth\talt_reads\t" + "\t".join(f"{r}_sfs\t{r}_cluster" for r in RUNS) + "\n")
        rows = []
        for c, p, ref, alt, t, dp, na in res:
            flags = []
            for r in RUNS:
                s, e = p - PAD, p + len(ref) + PAD
                flags += [int(hit(cov[r][0], c, s, e)), int(hit(cov[r][1], c, s, e))]
            rows.append((t, dp, na, flags))
            fo.write(f"{c}\t{p + 1}\t{ref}\t{alt}\t{t}\t{dp}\t{na}\t" + "\t".join(map(str, flags)) + "\n")

    print(f"# PAD={PAD} WIN={WIN} MIN_ALT={MIN_ALT}; variants inside all.bed: {len(rows)}")
    print("type\tn\tdepth_median\talt_median\tin_reads(alt>=MIN_ALT)\t"
          + "\t".join(f"{r}_sfs_all\t{r}_cluster_all\t{r}_sfs_inreads\t{r}_cluster_inreads" for r in RUNS))
    for t in ("SNV", "INDEL_HP", "INDEL_other", "MNV"):
        sel = [x for x in rows if x[0] == t]
        if not sel:
            continue
        n = len(sel)
        dps = sorted(x[1] for x in sel)
        nas = sorted(x[2] for x in sel)
        inr = [x for x in sel if x[2] >= MIN_ALT]
        cells = [t, str(n), str(dps[n // 2]), str(nas[n // 2]), f"{len(inr)} ({len(inr) / n:.3f})"]
        for k in range(len(RUNS)):
            a = sum(x[3][2 * k] for x in sel)
            b = sum(x[3][2 * k + 1] for x in sel)
            ai = sum(x[3][2 * k] for x in inr)
            bi = sum(x[3][2 * k + 1] for x in inr)
            cells += [f"{a / n:.3f}", f"{b / n:.3f}",
                      f"{ai / len(inr):.3f}" if inr else "-", f"{bi / len(inr):.3f}" if inr else "-"]
        print("\t".join(cells))
