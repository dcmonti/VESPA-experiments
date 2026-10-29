#!/usr/bin/env python3
"""Allele lengths, read by read, at the tandem-repeat length changes of the V0.5 ground truth (EVENTTYPE=CNV:TR) in the
clonal and subclonal regions, the 7 whose breakpoints no SFS of the full germline set covers (lost_points.tsv) and, for
comparison, the others; and at the other records with a breakpoint without SFS (two subclonal deletions).

The benchmark gives the somatic event as the size change between the germline and the somatic allele of the same
haplotype (README of V0.5: a normal with a 42 bp insertion and a tumor with a 28 bp deletion relative to GRCh38 give a
70 bp somatic deletion). So the expected picture is: on one haplotype, tumor allele = normal allele + SVLEN.

Locus: the TRF intervals (Severus annotation) that overlap the record extended by its confidence intervals (CIPOS,
CIEND) and by PAD bp, merged, plus FLANK bp on each side; without such an interval, the record +- WIN bp.
Allele of a read: primary alignment with MAPQ >= MIN_MAPQ that spans the whole locus; delta = read bases between the
two ends of the locus - locus length (0 = the GRCh38 length; > 0 longer; < 0 shorter). HP tag of the read (0 untagged).
BAMs: tumor and normal, raw and error-corrected (the reads the SFS come from are the corrected ones).

Output: tr_alleles.tsv (one row per read), tr_alleles_summary.tsv (per record, sample, read set and HP: number of
reads, median delta, and the delta values with their counts, most frequent first).
Usage: python tr_alleles.py [MIN_MAPQ=20]
"""
import csv, os, re, statistics, sys
from collections import Counter, defaultdict
import pysam

MIN_MAPQ = int(sys.argv[1]) if len(sys.argv) > 1 else 20
PAD, FLANK, WIN = 50, 20, 300
A = os.environ["HG008"]
R = f"{A}/reads_ec_align"
O = os.environ["HG008"] + "/reads_ec_align/sfs_signal"; os.makedirs(O, exist_ok=True)
GT = f"{A}/ground_truths/v0.5/GRCh38_HG008-T-V0.5_somatic-stvar_ALL.draftbenchmark.vcf"
TRF = f"{A}/refs/sev_extras/human_GRCh38_no_alt_analysis_set.trf.bed"
BAMS = [("tumor", "ec", f"{R}/align/tumor.ec.hp.bam"), ("normal", "ec", f"{R}/align/normal.ec.hp.bam"),
        ("tumor", "raw", f"{R}/align/tumor.raw.hp.bam"), ("normal", "raw", f"{R}/align/normal.raw.hp.bam")]


def info(s, k):
    m = re.search(rf"(?:^|;){k}=([^;]+)", s)
    return m.group(1) if m else ""


def qpos_at(aln, ref):
    """Query position (hard clips excluded, as query_sequence) aligned to reference position ref; inside a deletion, the
    query position of the next aligned base. None if ref is outside the aligned part."""
    if ref < aln.reference_start or ref >= aln.reference_end: return None
    rp, qp = aln.reference_start, 0
    for op, n in aln.cigartuples:
        if op in (0, 7, 8):
            if rp <= ref < rp + n: return qp + (ref - rp)
            rp += n; qp += n
        elif op in (2, 3):
            if rp <= ref < rp + n: return qp
            rp += n
        elif op in (1, 4): qp += n
    return None


def loci():
    """Records as above, as (id, chrom, locus start, locus end, svtype, svlen, lost, pos, end)."""
    rows = list(csv.DictReader(open(f"{O}/bp_points.tsv"), delimiter="\t"))
    # all CNV:TR records, plus every other record with a breakpoint that no SFS covers (the subclonal deletions)
    lost = {r["id"] for r in rows if r["in_cs"] == "1" and int(r["full_sfs"]) == 0}
    ids = {r["id"] for r in rows if r["eventtype"] == "CNV:TR" and r["in_cs"] == "1"} | lost
    trf = defaultdict(list)
    for l in open(TRF):
        f = l.split("\t"); trf[f[0]].append((int(f[1]), int(f[2])))
    out = []
    for l in open(GT):
        if l.startswith("#"): continue
        c = l.rstrip("\n").split("\t")
        if c[2] not in ids: continue
        ch, pos = c[0], int(c[1]) - 1
        end = int(info(c[7], "END")) - 1 if info(c[7], "END") else pos
        cp = [int(x) for x in info(c[7], "CIPOS").split(",")] if info(c[7], "CIPOS") else [0, 0]
        ce = [int(x) for x in info(c[7], "CIEND").split(",")] if info(c[7], "CIEND") else cp
        lo, hi = pos + cp[0] - PAD, end + ce[1] + PAD
        iv = sorted((s, e) for s, e in trf[ch] if s < hi and e > lo)
        if iv: s, e = min(x[0] for x in iv) - FLANK, max(x[1] for x in iv) + FLANK
        else: s, e = pos - WIN, end + WIN
        out.append((c[2], ch, s, e, info(c[7], "SVTYPE"), int(info(c[7], "SVLEN")), c[2] in lost, pos + 1, end + 1))
    return out


if __name__ == "__main__":
    L = loci()
    per = defaultdict(list)
    with open(f"{O}/tr_alleles.tsv", "w") as fo:
        fo.write("id\tchrom\tlocus_start\tlocus_end\tsample\treads\tread\tmapq\thp\tdelta\n")
        for sample, kind, path in BAMS:
            bam = pysam.AlignmentFile(path)
            for id_, ch, s, e, t, svlen, lost, pos, end in L:
                for a in bam.fetch(ch, s, e):
                    if a.is_secondary or a.is_supplementary or a.mapping_quality < MIN_MAPQ: continue
                    if a.reference_start > s or a.reference_end <= e: continue
                    # read bases from reference position s (included) to e (excluded): an insertion or deletion at either
                    # end of the locus falls in the flank, which is unique sequence
                    q1, q2 = qpos_at(a, s), qpos_at(a, e)
                    if q1 is None or q2 is None: continue
                    d = (q2 - q1) - (e - s)
                    hp = a.get_tag("HP") if a.has_tag("HP") else 0
                    per[(id_, sample, kind, hp)].append(d)
                    fo.write(f"{id_}\t{ch}\t{s}\t{e}\t{sample}\t{kind}\t{a.query_name}\t{a.mapping_quality}\t{hp}\t{d}\n")
    with open(f"{O}/tr_alleles_summary.tsv", "w") as fo:
        fo.write(f"# MIN_MAPQ={MIN_MAPQ} PAD={PAD} FLANK={FLANK}\n")
        fo.write("id\tlost\tsvtype\tsvlen\tchrom\tpos\tend\tlocus_start\tlocus_end\tsample\treads\thp\tn\tmedian_delta\tdeltas\n")
        for id_, ch, s, e, t, svlen, lost, pos, end in L:
            for sample, kind, _ in BAMS:
                for hp in (1, 2, 0):
                    v = per.get((id_, sample, kind, hp), [])
                    if not v: continue
                    hist = ",".join(f"{d}:{n}" for d, n in Counter(v).most_common())
                    fo.write(f"{id_}\t{int(lost)}\t{t}\t{svlen}\t{ch}\t{pos}\t{end}\t{s}\t{e}\t{sample}\t{kind}\t{hp}\t{len(v)}\t"
                             f"{statistics.median(v):g}\t{hist}\n")
    print("done", file=sys.stderr)
