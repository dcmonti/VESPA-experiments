#!/usr/bin/env python3
"""Per-breakpoint coverage of the V0.5 ground truth, for each index of the germline cascade (tab:res-cascade).

Same points, tolerance and stage definitions as bp_coverage/bp_point_table_prov.py (TOL 50), but one row per
breakpoint instead of the aggregate, so that the coverage can be split by SV class and calling path, and the
breakpoints lost between two indexes can be listed.

Indexes (columns of tab:res-cascade):
  ref   germline_cascade/rows/R1_ref_nofilters      SFS only (the call is not computable)
  utg   germline_cascade/rows/R2_refutg_nofilters
  full  bp_coverage (sfs_clean.bed.gz) + bp_coverage/nogl_w5/provenance.tsv
  v63   germline_cascade/rows/R2w_refv63_nofilters
Stages per index: sfs (placed SFS), clus (EVAL pass = 1), cand (sources of the CALL rows = candidate SVs),
split into cand_placed (S and F tokens: placed SFS, templated-insertion fragments) and cand_clip (C and J tokens:
clip clusters, deletions through, BND rescue, junctions between supplementary alignments).
For each SFS hit the read name is kept (sfs_hits.tsv.gz), for the diagnosis of the lost breakpoints.
A background rate is given by NRAND random points inside the clonal and subclonal regions (seed 20261014).

Output (in this directory): bp_points.tsv, sfs_hits.tsv.gz, bp_background.tsv
Usage: python bp_points.py
"""
import bisect, gzip, os, random, re, subprocess, sys
from collections import defaultdict

TOL = 50
NRAND = 20000
A = os.environ["HG008"]
R = f"{A}/reads_ec_align"
B = f"{A}/ground_truths/v0.5"
GT = f"{B}/GRCh38_HG008-T-V0.5_somatic-stvar_ALL.draftbenchmark.vcf"
TRF = f"{A}/refs/sev_extras/human_GRCh38_no_alt_analysis_set.trf.bed"
O = os.environ["HG008"] + "/reads_ec_align/sfs_signal"; os.makedirs(O, exist_ok=True)
IDX = [("ref", f"{R}/germline_cascade/rows/R1_ref_nofilters/sfs_clean.bed.gz", None),
       ("utg", f"{R}/germline_cascade/rows/R2_refutg_nofilters/sfs_clean.bed.gz",
        f"{R}/germline_cascade/rows/R2_refutg_nofilters/provenance.tsv"),
       ("full", f"{R}/bp_coverage/sfs_clean.bed.gz", f"{R}/bp_coverage/nogl_w5/provenance.tsv"),
       ("v63", f"{R}/germline_cascade/rows/R2w_refv63_nofilters/sfs_clean.bed.gz",
        f"{R}/germline_cascade/rows/R2w_refv63_nofilters/provenance.tsv")]


def load_bed(p):
    d = defaultdict(list)
    for l in open(p):
        c = l.split("\t"); d[c[0]].append((int(c[1]), int(c[2])))
    for k in d: d[k].sort()
    return d


def inside(reg, ch, p):
    v = reg.get(ch)
    if not v: return False
    j = bisect.bisect_right([x[0] for x in v], p) - 1
    return j >= 0 and v[j][1] > p


def near(d, ch, p, pad):
    v = d.get(ch)
    if not v: return False
    j = bisect.bisect_right([x[0] for x in v], p + pad)
    k = j - 1
    while k >= 0 and k > j - 2000:
        if v[k][1] > p - pad: return True
        k -= 1
    return False


def info(s, k):
    m = re.search(rf"(?:^|;){k}=([^;]+)", s)
    return m.group(1) if m else ""


def gt_points():
    """As bp_point_table_prov.gt_points, with the record attached: one point per BND record, POS of an INS,
    POS and END of DEL/DUP/INV; records shorter than 50 bp are skipped."""
    pts = []
    for l in open(GT):
        if l.startswith("#"): continue
        c = l.rstrip("\n").split("\t"); ch, pos = c[0], int(c[1]) - 1
        t = info(c[7], "SVTYPE")
        meta = dict(id=c[2], svtype=t, eventtype=info(c[7], "EVENTTYPE"), filter=c[6],
                    subclonal=info(c[7], "SUBCLONAL"), complex=info(c[7], "COMPLEX"),
                    vaf=info(c[7], "SVVIZ_VAF_ALL"), svlen="", mate="")
        if "[" in c[4] or "]" in c[4] or t == "BND":
            m = re.search(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]", c[4])
            meta["mate"] = f"{m.group(1)}:{m.group(2)}" if m else ""
            pts.append((ch, pos, "bnd", meta)); continue
        e = info(c[7], "END"); sl = info(c[7], "SVLEN")
        L = abs(int(sl)) if sl else (int(e) - pos if e else 0)
        if L < 50: continue
        meta["svlen"] = str(L)
        pts.append((ch, pos, "pos", meta))
        if t != "INS" and e: pts.append((ch, int(e) - 1, "end", meta))
    seen, out = set(), []
    for p in sorted(pts, key=lambda x: (x[0], x[1])):
        if (p[0], p[1]) in seen: continue
        seen.add((p[0], p[1])); out.append(p)
    return out


TOK = re.compile(r"(S|C|J|F):([^:|>+\t]+):(\d+)(?:-(\d+))?(?::[LR])?(?:>([^:|>+\t]+):(\d+))?")


def objs_of(s, with_frag):
    out = []
    for k, ch, a, b, mch, mp in TOK.findall(s):
        if k == "S": out.append((k, ((ch, int(a), int(b) + 1),)))
        elif k in "CJ":
            o = [(ch, int(a) - 1, int(a))]
            if mch: o.append((mch, int(mp) - 1, int(mp)))
            out.append((k, tuple(o)))
        elif k == "F" and with_frag:
            out.append((k, ((ch, int(a) - 1, int(a)), (ch, int(b) - 1, int(b)))))
    return out


def prov_sets(prov):
    """clus, cand, cand_placed, cand_clip as {chrom: sorted intervals}, as bp_point_table_prov.py."""
    sets = {k: defaultdict(list) for k in ("clus", "cand", "cand_placed", "cand_clip")}
    def add(name, o):
        for ch, s, e in o: sets[name][ch].append((s, e))
    for l in open(prov):
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        if f[0] == "EVAL" and f[-1] == "1":
            for k, o in objs_of(f[2], False): add("clus", o)
        elif f[0] == "CALL":
            for k, o in objs_of(f[8], True):
                add("cand", o); add("cand_placed" if k in "SF" else "cand_clip", o)
            for k, o in objs_of(re.sub(r"(^|\|)[SCJ]:[^|]*", "", f[8]), True): add("clus", o)
    for d in sets.values():
        for k in d: d[k].sort()
    return sets


def sfs_hits(path, wins):
    """Stream the SFS file; return {point index: [(start, end, read), ...]} for SFS within TOL of a point."""
    bins = defaultdict(list)
    for i, (ch, p) in enumerate(wins):
        for b in range((p - TOL) // 10000, (p + TOL) // 10000 + 1): bins[(ch, b)].append((i, p))
    hits = defaultdict(list)
    zc = subprocess.Popen(["zcat", path], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for l in zc.stdout:
        c = l.split("\t", 4); ch = c[0]; s = int(c[1]); e = int(c[2])
        for b in range(s // 10000, (e - 1) // 10000 + 1):
            for i, p in bins.get((ch, b), ()):
                if s < p + TOL + 1 and e > p - TOL: hits[i].append((s, e, c[3]))
    zc.wait()
    return hits


if __name__ == "__main__":
    regs = {r: load_bed(f"{B}/GRCh38_HG008-T-V0.5_somatic-stvar-{r}.draftbenchmark.bed")
            for r in ("clonal", "clonal_and_subclonal")}
    trf = load_bed(TRF)
    pts = gt_points()
    # random background points inside the clonal and subclonal regions
    rnd = random.Random(20261014); cs = regs["clonal_and_subclonal"]
    segs = [(ch, s, e) for ch, v in cs.items() for s, e in v]; cum = [0]
    for ch, s, e in segs: cum.append(cum[-1] + e - s)
    bg = []
    for _ in range(NRAND):
        x = rnd.randrange(cum[-1]); j = bisect.bisect_right(cum, x) - 1
        bg.append((segs[j][0], segs[j][1] + x - cum[j]))
    wins = [(ch, p) for ch, p, _, _ in pts] + bg
    res = {}; bgres = {}
    with gzip.open(f"{O}/sfs_hits.tsv.gz", "wt") as fh:
        fh.write("index\tpoint\tchrom\tpos\tsfs_start\tsfs_end\tread\n")
        for name, sfs, prov in IDX:
            print(f"[bp_points] {name}: SFS", file=sys.stderr, flush=True)
            h = sfs_hits(sfs, wins)
            sets = prov_sets(prov) if prov else None
            for i, (ch, p, role, m) in enumerate(pts):
                r = {"sfs": len(h.get(i, []))}
                for k in ("clus", "cand", "cand_placed", "cand_clip"):
                    r[k] = int(near(sets[k], ch, p, TOL)) if sets else ""
                res[(name, i)] = r
                for s, e, rd in h.get(i, []): fh.write(f"{name}\t{m['id']}:{role}\t{ch}\t{p + 1}\t{s}\t{e}\t{rd}\n")
            bgres[name] = sum(1 for j in range(len(pts), len(wins)) if h.get(j))
    with open(f"{O}/bp_points.tsv", "w") as fo:
        cols = ["id", "role", "chrom", "pos", "svtype", "eventtype", "svlen", "mate", "filter", "subclonal", "complex",
                "vaf", "in_clonal", "in_cs", "in_trf"]
        for name, _, _ in IDX:
            cols += [f"{name}_{k}" for k in ("sfs", "clus", "cand", "cand_placed", "cand_clip")]
        fo.write("\t".join(cols) + "\n")
        for i, (ch, p, role, m) in enumerate(pts):
            row = [m["id"], role, ch, str(p + 1), m["svtype"], m["eventtype"], m["svlen"], m["mate"], m["filter"],
                   m["subclonal"], m["complex"], m["vaf"], str(int(inside(regs["clonal"], ch, p))),
                   str(int(inside(regs["clonal_and_subclonal"], ch, p))), str(int(near(trf, ch, p, TOL)))]
            for name, _, _ in IDX:
                r = res[(name, i)]
                row += [str(r[k]) for k in ("sfs", "clus", "cand", "cand_placed", "cand_clip")]
            fo.write("\t".join(row) + "\n")
    with open(f"{O}/bp_background.tsv", "w") as fo:
        fo.write("index\trandom_points\twith_sfs_within_50bp\tfraction\n")
        for name, _, _ in IDX:
            fo.write(f"{name}\t{NRAND}\t{bgres[name]}\t{bgres[name] / NRAND:.4f}\n")
    print("done", file=sys.stderr)
