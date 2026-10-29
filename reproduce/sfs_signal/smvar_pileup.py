#!/usr/bin/env python3
"""Alternative reads of the clonal somatic small variants (GIAB smvar v0.3) in the tumor and normal EC alignments
of reads_ec_align. Redo of het6/hpc_caps/sfs_strat/{small_variant_support,normal_check}.py (2026-09-25).

Variants: those of smvar_cov.tsv.gz (sfs_classes.py: inside smvar all.bed x V0.5 clonal and subclonal regions).
Counting, unchanged: primary alignments with MAPQ >= 60. SNV: the read base at POS equals ALT. Indel: the read has
an indel of the same type and length starting within +-WIN bp of POS (tolerates left/right placement in repeats).
A variant is "in our reads" with at least MIN_ALT tumor alternative reads.

Output: smvar_pileup.tsv.gz (per variant: tumor and normal depth and alternative reads), smvar_pileup_summary.tsv
Usage: python smvar_pileup.py [THREADS=24] [CHROMS=all, comma-separated, for a test]
"""
import gzip, os, sys
from collections import defaultdict
from multiprocessing import Pool
import pysam

THREADS = int(sys.argv[1]) if len(sys.argv) > 1 else 24
CHROMS = set(sys.argv[2].split(",")) if len(sys.argv) > 2 else None
WIN, MIN_ALT = 5, 3
R = os.environ["HG008"] + "/reads_ec_align"
BAMS = {"tumor": f"{R}/align/tumor.ec.hp.bam", "normal": f"{R}/align/normal.ec.hp.bam"}
O = os.environ["HG008"] + "/reads_ec_align/sfs_signal"; os.makedirs(O, exist_ok=True)


def count_alt(bam, c, p, ref, alt, t):
    dp = na = 0
    if t == "SNV" and len(ref) == 1:
        for col in bam.pileup(c, p, p + 1, truncate=True, min_mapping_quality=60, stepper="nofilter",
                              ignore_orphans=False, min_base_quality=0):
            for pr in col.pileups:
                a = pr.alignment
                # pileup(min_mapping_quality=60) does not drop these reads with stepper="nofilter" (pysam 0.24.0): filter here
                if a.is_secondary or a.is_supplementary or a.mapping_quality < 60: continue
                dp += 1
                if not pr.is_del and not pr.is_refskip and a.query_sequence[pr.query_position].upper() == alt.upper():
                    na += 1
        return dp, na
    if t == "SNV": return -1, -1  # MNV: not counted
    ins = len(alt) > len(ref); L = abs(len(alt) - len(ref))
    for a in bam.fetch(c, p, p + 1):
        if a.is_secondary or a.is_supplementary or a.mapping_quality < 60: continue
        if a.reference_start > p - WIN or a.reference_end < p + len(ref) + WIN: continue
        dp += 1; rp = a.reference_start
        for op, n in a.cigartuples:
            if op in (0, 7, 8): rp += n
            elif op == 1:
                if ins and n == L and abs(rp - (p + 1)) <= WIN: na += 1; break
            elif op == 2:
                if not ins and n == L and abs(rp - (p + 1)) <= WIN: na += 1; break
                rp += n
            elif op == 3: rp += n
            if rp > p + WIN + len(ref): break
    return dp, na


def work(job):
    c, vs = job; out = []
    bams = {k: pysam.AlignmentFile(v) for k, v in BAMS.items()}
    for p, ref, alt, t in vs:
        out.append((c, p, ref, alt, t) + count_alt(bams["tumor"], c, p, ref, alt, t)
                   + count_alt(bams["normal"], c, p, ref, alt, t))
    return out


if __name__ == "__main__":
    byc = defaultdict(list)
    with gzip.open(f"{O}/smvar_cov.tsv.gz", "rt") as fh:
        next(fh)
        for l in fh:
            f = l.split("\t")
            if CHROMS and f[0] not in CHROMS: continue
            byc[f[0]].append((int(f[1]) - 1, f[2], f[3], f[4]))
    jobs = [(c, vs[i:i + 200]) for c, vs in byc.items() for i in range(0, len(vs), 200)]
    with Pool(THREADS) as pool:
        res = [x for r in pool.imap_unordered(work, jobs) for x in r]
    suf = "" if CHROMS is None else ".test"
    with gzip.open(f"{O}/smvar_pileup{suf}.tsv.gz", "wt") as fo:
        fo.write("chrom\tpos\tref\talt\ttype\ttumor_depth\ttumor_alt\tnormal_depth\tnormal_alt\n")
        for r in sorted(res): fo.write(f"{r[0]}\t{r[1] + 1}\t" + "\t".join(map(str, r[2:])) + "\n")
    with open(f"{O}/smvar_pileup_summary{suf}.tsv", "w") as fo:
        fo.write("type\tn\ttumor_depth_med\ttumor_alt_med\tin_reads\tn_normal_dp10\tnormal_alt0\t"
                 "nvaf_0_5\tnvaf_5_20\tnvaf_20_35\tnvaf_35_65\tnvaf_65_100\ttumor_vaf_med\n")
        for t in ("SNV", "INDEL_HP", "INDEL_other"):
            sel = [r for r in res if r[4] == t and r[5] >= 0]
            if not sel: continue
            n = len(sel); med = lambda v: sorted(v)[len(v) // 2]
            inr = sum(1 for r in sel if r[6] >= MIN_ALT)
            nn = [r for r in sel if r[7] >= 10]; m = len(nn)
            vaf = [r[8] / r[7] for r in nn]
            bins = [sum(1 for v in vaf if lo <= v < hi) / m for lo, hi in ((0, .05), (.05, .2), (.2, .35), (.35, .65), (.65, 1.01))]
            tv = med([r[6] / r[5] for r in sel if r[5] > 0])
            fo.write(f"{t}\t{n}\t{med([r[5] for r in sel])}\t{med([r[6] for r in sel])}\t{inr / n:.4f}\t{m}\t"
                     f"{sum(1 for r in nn if r[8] == 0) / m:.4f}\t" + "\t".join(f"{b:.4f}" for b in bins) + f"\t{tv:.3f}\n")
