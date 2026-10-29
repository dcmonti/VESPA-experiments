#!/usr/bin/env python3
"""Truth coverage with per-read attribution (possible only on simulated data: the name of every read
tells which sequence it comes from). One cell; writes cells/<S>/<EC>/eval/attr.tsv and attr_background.tsv.
  Short variants: covered at the SFS stage if a projected SFS (prov/sfs_clean.bed.gz) contains the position (SNV:
    the base; indel: [POS-5, POS+len(REF)+5], for the shift in repeats) and comes from a read that carries
    the alternative allele in the corrected BAM (SNV: base = ALT; indel: same type and length within +-5 bp;
    primary alignments, MAPQ 60, as in small_variant_support.py).
  SV breakpoints (2026-09-29, rule b): covered at the SFS stage if an SFS within 50 bp comes from a
    read that (1) comes from a tumor sequence carrying the SV (hap1_T1 or hap1_T2 for the clonal ones, only hap1_T2
    for the subclonal ones) and (2) carries the SV at that point: one of its alignments, primary or supplementary, has
    a clipped end within 50 bp of the breakpoint, or an indel of at least MIN_OP bp in the CIGAR within
    50 bp. The projected SFS include primary and supplementary alignments (sfs_to_genome.py): a read with a
    junction is often split into several alignments.
    Tandem DUP (2026-09-29): a read that contains both copies is aligned with
    an insertion of the copy at any point of the duplicated segment, and the junction between the copies (and its
    SFS) ends up there. For both breakpoints of a DUP the window is therefore the whole segment +-50 bp
    (signal, SFS and cluster), instead of +-50 bp from the breakpoint.
  Cluster: covered if covered at the SFS stage and a cluster of weight >= 5 (EVAL pass=1 of the dump at m_min 5) lies within the
    same tolerance (10 bp for short variants, 50 bp for breakpoints).
  Background: fraction of 5 000 random positions (>= 2 kb from every variant) that are covered under the same
    rule as the breakpoints (read of hap1_T1 or hap1_T2 with clip or indel >= MIN_OP within 50 bp and an SFS within 50 bp),
    and, for comparison, under condition (1) alone.
USAGE: python3 07_attrib.py D S EC"""
import sys, os, re, gzip, bisect, random
from collections import defaultdict
import pysam
D, S, EC = sys.argv[1:4]
W = f"{D}/cells/{S}/{EC}"; E = f"{W}/eval"; os.makedirs(E, exist_ok=True)
CTG = "sim"; WIN = 5; PAD_SMALL = 10; TOL_BP = 50; MIN_OP = 20
comp = lambda r: r.split(".")[-2]                       # <sample>.<component>.<n>
TUMOR = {"clonal": {"hap1_T1", "hap1_T2"}, "subclonal": {"hap1_T2"}}
KEY = {"truncal": "clonal", "subclonal": "subclonal"}   # internal key of the truth files

# --- projected SFS, per read and as intervals ---------------------------------------------------
by_read = defaultdict(list); tum_iv = []
for l in gzip.open(f"{W}/prov/sfs_clean.bed.gz", "rt"):
    f = l.split("\t")
    if f[0] != CTG: continue
    s, e, r = int(f[1]), int(f[2]), f[3]
    by_read[r].append((s, e))
    if comp(r) in TUMOR["clonal"]: tum_iv.append((s, e, comp(r)))
tum_iv.sort(); tum_start = [x[0] for x in tum_iv]; tum_maxlen = max((e - s for s, e, _ in tum_iv), default=0)
def sv_reads(bam, p, comps, lo=None, hi=None):
    """reads of the sequences comps with an SV signal within TOL_BP of p: clipped end (soft/hard) of a
    primary or supplementary alignment, or an indel >= MIN_OP in the CIGAR"""
    lo = p - TOL_BP if lo is None else lo; hi = p + TOL_BP if hi is None else hi   # window [lo, hi]
    out = set()
    for a in bam.fetch(CTG, max(0, lo - 1), hi + 1):
        if a.is_secondary or a.is_unmapped or comp(a.query_name) not in comps: continue
        cg = a.cigartuples
        if cg[0][0] in (4, 5) and lo <= a.reference_start <= hi: out.add(a.query_name); continue
        if cg[-1][0] in (4, 5) and lo <= a.reference_end <= hi: out.add(a.query_name); continue
        rp = a.reference_start
        for op, n in cg:
            if op in (1, 2) and n >= MIN_OP and rp <= hi and rp + (n if op == 2 else 0) >= lo:
                out.add(a.query_name); break
            if op in (0, 2, 3, 7, 8): rp += n
            if rp > hi: break
    return out
def attributed_near(p, reads, lo=None, hi=None):
    lo = p - TOL_BP if lo is None else lo; hi = p + TOL_BP if hi is None else hi
    return any(a < hi + 1 and b > lo for r in reads for a, b in by_read.get(r, ()))
def tumor_sfs_near(p, comps):
    i = bisect.bisect_right(tum_start, p + TOL_BP); k = i - 1
    while k >= 0 and tum_iv[k][0] > p - TOL_BP - tum_maxlen:
        s, e, c = tum_iv[k]
        if e > p - TOL_BP and c in comps: return True
        k -= 1
    return False
cl = []
TOK = re.compile(r"(S|C|J):([^:|>+\t]+):(\d+)(?:-(\d+))?(?::[LR])?(?:>([^:|>+\t]+):(\d+))?")
for l in open(f"{W}/prov/provenance.tsv"):
    f = l.rstrip("\n").split("\t")
    if f[0] != "EVAL" or f[-1] != "1": continue
    for k, ch, a, b, mch, mp in TOK.findall(f[2]):
        if k == "S" and ch == CTG: cl.append((int(a), int(b) + 1))
        elif k in "CJ":
            if ch == CTG: cl.append((int(a) - 1, int(a)))
            if mch == CTG: cl.append((int(mp) - 1, int(mp)))
cl.sort(); cl_start = [a for a, _ in cl]; cl_maxlen = max((b - a for a, b in cl), default=0)
def cluster_near(s, e):
    i = bisect.bisect_right(cl_start, e); k = i - 1
    while k >= 0 and cl[k][0] > s - cl_maxlen:
        if cl[k][1] > s: return True
        k -= 1
    return False

# --- short variants --------------------------------------------------------------------------------
def alt_reads(bam, p, ref, alt):
    """names of the reads (primary, MAPQ 60) with the alternative allele"""
    out = set()
    if len(ref) == len(alt) == 1:
        for col in bam.pileup(CTG, p, p + 1, truncate=True, min_mapping_quality=60, stepper="nofilter", min_base_quality=0):
            for pr in col.pileups:
                a = pr.alignment
                if a.is_secondary or a.is_supplementary or pr.is_del or pr.is_refskip: continue
                if a.query_sequence[pr.query_position].upper() == alt.upper(): out.add(a.query_name)
        return out
    ins = len(alt) > len(ref); L = abs(len(alt) - len(ref))
    for a in bam.fetch(CTG, max(0, p - WIN), p + len(ref) + WIN):
        if a.is_secondary or a.is_supplementary or a.mapping_quality < 60: continue
        rp = a.reference_start
        for op, n in a.cigartuples:
            if op in (0, 7, 8): rp += n
            elif op == 1:
                if ins and n == L and abs(rp - (p + 1)) <= WIN: out.add(a.query_name); break
            elif op == 2:
                if not ins and n == L and abs(rp - (p + 1)) <= WIN: out.add(a.query_name); break
                rp += n
            elif op == 3: rp += n
            if rp > p + WIN + len(ref): break
    return out
res = defaultdict(lambda: [0, 0, 0])      # (type, clone) -> [n, SFS, cluster]
with pysam.AlignmentFile(f"{W}/tumor.ec.bam") as bam:
    for l in open(f"{D}/truth/truth_small.vcf"):
        if l[0] == "#": continue
        f = l.rstrip("\n").split("\t"); p = int(f[1]) - 1; ref, alt = f[3], f[4]
        info = dict(x.split("=") for x in f[7].split(";") if "=" in x)
        key = (info["VTYPE"], KEY[info["CLONE"]]); res[key][0] += 1
        snv = len(ref) == len(alt) == 1
        s, e = (p, p + 1) if snv else (p - WIN, p + len(ref) + WIN)
        ok = any(a < e and b > s for r in alt_reads(bam, p, ref, alt) for a, b in by_read.get(r, ()))
        if ok:
            res[key][1] += 1
            if cluster_near(p - PAD_SMALL, p + len(ref) + PAD_SMALL): res[key][2] += 1

# --- SV breakpoints ------------------------------------------------------------------------------------
bam = pysam.AlignmentFile(f"{W}/tumor.ec.bam")
for l in open(f"{D}/truth/truth_sv.vcf"):
    if l[0] == "#": continue
    f = l.rstrip("\n").split("\t"); pos = int(f[1]) - 1
    info = dict(x.split("=") for x in f[7].split(";") if "=" in x)
    c = KEY[info["CLONE"]]; pts = [pos] if info["SVTYPE"] == "INS" else [pos, int(info["END"]) - 1]
    # DUP: window = duplicated segment +-50 bp, for both breakpoints
    lo, hi = (pos - TOL_BP, int(info["END"]) + TOL_BP) if info["SVTYPE"] == "DUP" else (None, None)
    for p in pts:
        key = ("SV_bp", c); res[key][0] += 1
        if attributed_near(p, sv_reads(bam, p, TUMOR[c], lo, hi), lo, hi):
            res[key][1] += 1
            if cluster_near(p - TOL_BP if lo is None else lo, (p + TOL_BP if hi is None else hi) + 1): res[key][2] += 1

with open(f"{E}/attr.tsv", "w") as o:
    o.write("sample\tec\ttype\tclone\tn\tcovered_SFS\tcovered_cluster_w5\n")
    for ty in ("SNV", "INDEL_other", "SV_bp"):
        for c in ("clonal", "subclonal"):
            n, a, b = res[(ty, c)]
            if n: o.write(f"{S}\t{EC}\t{ty}\t{c}\t{n}\t{a}\t{b}\n")

# --- background ------------------------------------------------------------------------------------
truth = sorted([int(l.split("\t")[1]) - 1 for fn in ("truth_small.vcf", "truth_sv.vcf")
                for l in open(f"{D}/truth/{fn}") if l[0] != "#"])
L = int(open(f"{D}/ref/region.fa.fai").read().split("\t")[1]); rng = random.Random(1); pts = []
while len(pts) < 5000:
    p = rng.randrange(10000, L - 10000); i = bisect.bisect_left(truth, p)
    if (i < len(truth) and truth[i] - p < 2000) or (i > 0 and p - truth[i - 1] < 2000): continue
    pts.append(p)
bg = sum(attributed_near(p, sv_reads(bam, p, TUMOR["clonal"])) for p in pts) / len(pts)
bg1 = sum(tumor_sfs_near(p, TUMOR["clonal"]) for p in pts) / len(pts)
with open(f"{E}/attr_background.tsv", "w") as o:
    o.write("sample\tec\tbackground_rule_b\tbackground_tumor_SFS_only\n" f"{S}\t{EC}\t{bg:.4f}\t{bg1:.4f}\n")
print(f"{S} {EC}: attribution written, background rule b {bg:.4f} (tumor reads only {bg1:.4f})")
