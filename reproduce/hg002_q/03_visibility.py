#!/usr/bin/env python3
"""HG002, chr1: for every germline SV of the truth, which spanning reads of each homolog show it in their alignment.

Events: the DEL and INS of the dipcall VCF of T2T-HG002-Q100 v1.1 (GIAB v5.0q) with a length of 50 bp to MAX_LEN
(SVs that the reads traverse), near-pure indels (shorter allele at most 10 bp), genotype 0|1, 1|0 or 1|1, entirely
inside the benchmark regions of the SV set. Carrying homolog: 1 for 1|0, 2 for 0|1, both for 1|1.

Reads, as VESPA examines the normal (Caller::is_germline and the stats pass of caller.cpp): primary alignments with
MAPQ >= MIN_MAPQ that cover the breakpoint s (pos <= s <= end - 1). HP from 02_haplotag.sh (truth phasing).

A read SHOWS the event under two criteria, both copied from VESPA v1.0.1, on the alignment as it is (no ksw2
realignment of the normal read, which could only add reads):
  strict  the same-allele test of is_germline(): an op of the same type starting in [s - 150, s + len + 150] with
          length nl such that (nl >= 50 and min/max >= 0.90) or (nl >= 10 and |len - nl| < 50)
  loose   the length test of GERMLINE (gate A): an op of the same type, >= 50 bp, starting in [s - 500, s + len + 500],
          with |nl - len| <= max(50, 0.25 len)
(s + len bounds the window for INS too, as in the code.)

Output: per-SV TSV with, for each HP class (1, 2, 0 = untagged), the spanning reads and those that show the event
under each criterion. Usage: python 03_visibility.py OUT_DIR [CHROM [TAG]]   (TAG: the alignment, hg002.TAG.raw.hp.bam)
"""
import bisect, gzip, os, sys
from collections import defaultdict
import pysam

OUT = sys.argv[1]
CHROM = sys.argv[2] if len(sys.argv) > 2 else "chr1"
TAG = sys.argv[3] if len(sys.argv) > 3 else CHROM
MIN_LEN, MAX_LEN, MIN_MAPQ = 50, 10000, 60
STRICT_WIN, LOOSE_WIN, MIN_RATIO, DIFF_MAX, DIFF_MIN_LEN, LOOSE_FRAC = 150, 500, 0.90, 50, 10, 0.25
TRF = os.environ["REFS"] + "/human_GRCh38_no_alt_analysis_set.trf.bed"
DIP = f"{OUT}/truth/GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.vcf.gz"
BED = f"{OUT}/truth/HG002_GRCh38_v5.0q_stvar.benchmark.bed"
BAM = f"{OUT}/align/hg002.{TAG}.raw.hp.bam"
OUTD = OUT if TAG == CHROM else f"{OUT}/{TAG}"   # per-chromosome results of a multi-chromosome run
TR_FLANK = 50


def intervals(path, chrom):
    iv = []
    op = gzip.open if path.endswith(".gz") else open
    for l in op(path, "rt"):
        f = l.split("\t")
        if f[0] == chrom:
            iv.append((int(f[1]), int(f[2])))
    iv.sort()
    merged = []
    for s, e in iv:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [m[0] for m in merged], [m[1] for m in merged]


def inside(iv, s, e):
    starts, ends = iv
    i = bisect.bisect_right(starts, s) - 1
    return i >= 0 and ends[i] >= e


def overlaps(iv, s, e):
    starts, ends = iv
    i = bisect.bisect_right(starts, e) - 1
    return i >= 0 and ends[i] > s


def shows(aln, want, s, L):
    """(strict, loose) for one read, walking its CIGAR as is_germline() and the stats pass do."""
    strict = loose = False
    rpos = aln.reference_start
    for op, n in aln.cigartuples:
        if op == want:
            if s - STRICT_WIN <= rpos <= s + L + STRICT_WIN:
                ratio_ok = n >= MIN_LEN and min(n, L) / max(n, L) >= MIN_RATIO
                diff_ok = n >= DIFF_MIN_LEN and abs(L - n) < DIFF_MAX
                strict = strict or ratio_ok or diff_ok
            if s - LOOSE_WIN <= rpos <= s + L + LOOSE_WIN and n >= MIN_LEN:
                loose = loose or abs(n - L) <= max(MIN_LEN, LOOSE_FRAC * L)
        if op in (0, 2, 3, 7, 8):
            rpos += n
        if rpos > s + L + LOOSE_WIN:
            break
    return strict, loose


def main():
    os.makedirs(OUTD, exist_ok=True)
    bench, trf = intervals(BED, CHROM), intervals(TRF, CHROM)
    bam = pysam.AlignmentFile(BAM)
    out = open(f"{OUTD}/visibility.{CHROM}.tsv", "w")
    cols = ["chrom", "pos", "type", "len", "gt", "carrier", "in_trf"]
    for hp in ("hp1", "hp2", "hp0"):
        cols += [f"span_{hp}", f"strict_{hp}", f"loose_{hp}"]
    out.write("\t".join(cols) + "\n")
    skipped = defaultdict(int)
    for rec in pysam.VariantFile(DIP).fetch(CHROM):
        gt = rec.samples[0]["GT"]
        if not rec.samples[0].phased or gt not in ((0, 1), (1, 0), (1, 1)) or len(rec.alts) != 1:
            skipped["genotype"] += 1
            continue
        r, a = len(rec.ref), len(rec.alts[0])
        L = abs(r - a)
        if L < MIN_LEN or L > MAX_LEN:
            continue
        if min(r, a) > 10:
            skipped["complex"] += 1
            continue
        typ, want = ("DEL", 2) if r > a else ("INS", 1)
        s = rec.pos            # VCF POS, the base before the event (0-based start of the event)
        if not inside(bench, s - 1, s + (L if typ == "DEL" else 1)):
            skipped["outside benchmark regions"] += 1
            continue
        carrier = {(1, 0): "1", (0, 1): "2", (1, 1): "both"}[gt]
        n = {k: [0, 0, 0] for k in ("hp1", "hp2", "hp0")}
        for aln in bam.fetch(CHROM, s, s + 1):
            if aln.is_secondary or aln.is_supplementary or aln.is_unmapped or aln.mapping_quality < MIN_MAPQ:
                continue
            if aln.reference_start > s or aln.reference_end - 1 < s:
                continue
            hp = aln.get_tag("HP") if aln.has_tag("HP") else 0
            k = f"hp{hp}" if hp in (1, 2) else "hp0"
            st, lo = shows(aln, want, s, L)
            n[k][0] += 1; n[k][1] += st; n[k][2] += lo
        tr = int(overlaps(trf, s - TR_FLANK, s + (L if typ == "DEL" else 1) + TR_FLANK))
        row = [CHROM, s, typ, L, "|".join(map(str, gt)), carrier, tr] + [x for k in ("hp1", "hp2", "hp0") for x in n[k]]
        out.write("\t".join(map(str, row)) + "\n")
    out.close()
    print("skipped:", dict(skipped), file=sys.stderr)


if __name__ == "__main__":
    main()
