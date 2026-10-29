#!/usr/bin/env python3
"""Coverage table per GT breakpoint, with the "clusters" and "candidate clusters" rows from the
--call-provenance dump of VESPA instead of the --verbose log and the position-based
attribution. GT points, tolerance, n and genome fraction as in bp_point_table.py.

Stages:
  placed SFS          SFSBP_DIR/sfs_clean.bed.gz (unchanged);
  clusters            every EVAL object with pass = 1 of the run without the normal: SFS subclusters
                      (interval), clip clusters, rescued DEL/BND and junctions between supplementaries
                      (two points, p and j, for each clip cluster or junction composing them),
                      plus the fragments of the templated insertions of the CALL rows;
  candidate clusters  union of the sources of the CALL rows of the run without the normal (same
                      representation), plus the fragment (two points) of the templated insertions.
Provenance tokens: S:chrom:s-e (0-based, e included), C:chrom:p:L|R[>chrom:j] (1-based),
J:chrom:p:L|R>chrom:j (1-based), F:chrom:s-e (1-based).
PASS_ONLY=1: only the CALL rows with FILTER PASS ("in-call, post-call" row of the cascade); the clusters do not change.
Added columns (2026-10-08, cascade table): n_bed = objects with at least one point inside the regions;
area_bed_mb = union of the objects (no padding) clipped to the regions, in Mb; a clip cluster counts
as its two 1 bp points. n_sv_bed (candidate clusters stage only) = candidate SVs: CALL rows, excluding the
"mate" records (second record of a BND pair), with at least one source point inside the regions.
usage: bp_point_table_prov.py SFSBP_DIR PROV_TSV LABEL   (env GT, REGIONS, TOL as in bp_point_table.py)"""
import sys, os, re, bisect
from collections import defaultdict
S, PROV, LAB = sys.argv[1:4]
TOL = int(os.environ.get("TOL", 10))
PASS_ONLY = os.environ.get("PASS_ONLY", "0") == "1"
A = os.environ["HG008"]; B = f"{A}/ground_truths/v0.5"
GT = os.environ.get("GT", f"{B}/GRCh38_HG008-T-V0.5_somatic-stvar_ALL.draftbenchmark.vcf")
REGIONS = os.environ.get("REGIONS", "clonal,clonal_and_subclonal").split(",")
FAI = os.environ.get("FAI", f"{A}/refs/hg38.primary.fa.gz.fai")
G = sum(int(l.split("\t")[1]) for l in open(FAI))
import gzip
def load_bed(p):
    d = defaultdict(list)
    op = gzip.open if p.endswith(".gz") else open
    for l in op(p, "rt"):
        c = l.split("\t"); d[c[0]].append((int(c[1]), int(c[2])))
    for k in d: d[k].sort()
    return d
def prep_reg(reg):
    reg["_s"] = {k: [x[0] for x in v] for k, v in reg.items() if k != "_s"}; return reg
def inside(reg, ch, p):
    v = reg.get(ch)
    if not v: return False
    j = bisect.bisect_right(reg["_s"][ch], p) - 1
    return j >= 0 and v[j][1] > p
def gt_points():
    pts = []
    for l in open(GT):
        if l.startswith("#"): continue
        c = l.split("\t"); ch, pos = c[0], int(c[1]) - 1
        t = re.search(r"SVTYPE=(\w+)", c[7]).group(1)
        if "[" in c[4] or "]" in c[4] or t == "BND": pts.append((ch, pos)); continue
        e = re.search(r"(?:^|;)END=(\d+)", c[7]); sl = re.search(r"(?:^|;)SVLEN=(-?\d+)", c[7])
        L = abs(int(sl.group(1))) if sl else (int(e.group(1)) - pos if e else 0)
        if L < 50: continue
        pts.append((ch, pos))
        if t != "INS" and e: pts.append((ch, int(e.group(1)) - 1))
    return sorted(set(pts))
TOK = re.compile(r"(S|C|J|F):([^:|>+\t]+):(\d+)(?:-(\d+))?(?::[LR])?(?:>([^:|>+\t]+):(\d+))?")
def objs_of(s, with_frag):
    """Objects (as lists of intervals) named in a provenance string."""
    out = []
    for k, ch, a, b, mch, mp in TOK.findall(s):
        if k == "S": out.append(((ch, int(a), int(b) + 1),))
        elif k in "CJ":
            o = [(ch, int(a) - 1, int(a))]
            if mch: o.append((mch, int(mp) - 1, int(mp)))
            out.append(tuple(o))
        elif k == "F" and with_frag:
            out.append(((ch, int(a) - 1, int(a)), (ch, int(b) - 1, int(b))))
    return out
clus, cand = set(), set(); cand_sv = []
kinds = defaultdict(int); origins = defaultdict(int); nosrc = []
for l in open(PROV):
    if l.startswith("#"): continue
    f = l.rstrip("\n").split("\t")
    if f[0] == "EVAL" and f[-1] == "1":
        kinds[f[1]] += 1
        for o in objs_of(f[2], False): clus.add(o)
    elif f[0] == "CALL":
        if PASS_ONLY and f[7] != "PASS": continue
        origins[f[8].split("|")[0]] += 1
        os_ = objs_of(f[8], True)
        if f[8].split("|")[0] != "mate": cand_sv.append(os_)
        if not os_: nosrc.append(f[1])
        for o in os_: cand.add(o)
        # the fragment of a templated insertion is an object of its own (threshold on its own reads,
        # TEMPL_MIN_READS, not on --min-cluster-weight): it also enters the "clusters" row, so
        # that every candidate cluster is also a cluster
        for o in objs_of(re.sub(r"(^|\|)[SCJ]:[^|]*", "", f[8]), True): clus.add(o)
def to_d(objs):
    d = defaultdict(list)
    for o in objs:
        for ch, s, e in o: d[ch].append((s, e))
    for k in d: d[k].sort()
    d["_s"] = {k: [x[0] for x in v] for k, v in d.items() if not k.startswith("_")}
    d["_maxlen"] = {k: max(e - s for s, e in v) for k, v in d.items() if not k.startswith("_")}
    return d
def covered(objs, ch, p):
    v = objs.get(ch)
    if not v: return False
    j = bisect.bisect_right(objs["_s"][ch], p + TOL); k = j - 1
    while k >= 0 and v[k][0] > p - TOL - objs["_maxlen"][ch]:
        if v[k][1] > p - TOL: return True
        k -= 1
    return False
def genome_frac(d):
    tot = 0
    for k, v in d.items():
        if k.startswith("_"): continue
        cs = ce = None
        for s, e in v:
            if cs is None or s > ce:
                if cs is not None: tot += ce - cs
                cs, ce = s, e
            else: ce = max(ce, e)
        if cs is not None: tot += ce - cs
    return 100 * tot / G
def in_bed(d, reg):
    """objects with at least one point (start of a segment) inside the regions, and area of the union in the regions"""
    n = 0; segs = defaultdict(list)
    for o in d:
        if any(inside(reg, ch, s) for ch, s, e in o): n += 1
        for ch, s, e in o: segs[ch].append((s, e))
    tot = 0
    for ch, v in segs.items():
        R = reg.get(ch)
        if not R: continue
        v.sort(); m = []
        for s, e in v:
            if m and s <= m[-1][1]: m[-1][1] = max(m[-1][1], e)
            else: m.append([s, e])
        i = j = 0
        while i < len(m) and j < len(R):
            s, e = max(m[i][0], R[j][0]), min(m[i][1], R[j][1])
            if s < e: tot += e - s
            if m[i][1] < R[j][1]: i += 1
            else: j += 1
    return n, tot / 1e6
sfs = load_bed(f"{S}/sfs_clean.bed.gz")
sfs["_s"] = {k: [x[0] for x in v] for k, v in sfs.items() if not k.startswith("_")}
sfs["_maxlen"] = {k: max(e - s for s, e in v) for k, v in sfs.items() if not k.startswith("_")}
stages = [("placed SFS", sfs, sum(len(v) for k, v in sfs.items() if not k.startswith("_"))),
          ("clusters", to_d(clus), len(clus)), ("candidate clusters", to_d(cand), len(cand))]
objsets = {"placed SFS": [((ch, s, e),) for ch, v in sfs.items() if not ch.startswith("_") for s, e in v],
           "clusters": clus, "candidate clusters": cand}
print(f"# {LAB}: EVAL pass=1 by kind {dict(kinds)}; CALL by origin {dict(origins)}; "
      f"CALL without source {len(nosrc)} {nosrc[:5]}")
pts_all = gt_points()
print("label\tbed\tstage\tn_genome\tgenome_pct\tpoints\tpoints_cov\trecall\tn_bed\tarea_bed_mb\tn_sv_bed")
for bed in REGIONS:
    if bed == "-": bed = "genome"; pts = list(pts_all)
    else:
        reg = prep_reg(load_bed(f"{B}/GRCh38_HG008-T-V0.5_somatic-stvar-{bed}.draftbenchmark.bed"))
        pts = [p for p in pts_all if inside(reg, *p)]
    for name, d, n in stages:
        cov = sum(1 for ch, p in pts if covered(d, ch, p))
        nb, ab = in_bed(objsets[name], reg) if bed != "genome" else (n, float("nan"))
        nsv = (sum(1 for os_ in cand_sv if any(inside(reg, ch, s_) for o in os_ for ch, s_, e_ in o))
               if name == "candidate clusters" and bed != "genome" else "")
        print(f"{LAB}\t{bed}\t{name}\t{n}\t{genome_frac(d):.3f}\t{len(pts)}\t{cov}\t{cov/len(pts):.4f}\t{nb}\t{ab:.3f}\t{nsv}")
