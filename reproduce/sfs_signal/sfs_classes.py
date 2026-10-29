#!/usr/bin/env python3
"""Somatic classes of the placed SFS and of the SFS clusters, and SFS/cluster coverage of the somatic small variants,
for the indexes of the germline cascade. Redo, on reads_ec_align, of het6/hpc_caps/sfs_strat/{sfs_somatic_classes,
small_variant_support}.py (2026-09-25), with the inputs of tab:res-cascade.

Regions: the V0.5 clonal and subclonal SV regions, as in tab:res-cascade. The smvar v0.3 regions exclude the
neighbourhood of the SVs, so an object that is not on an SV and starts outside them has no small-variant truth and
gets its own class.
Objects:
  placed SFS   rows of sfs_clean.bed.gz (one row per placement, as in tab:res-cascade), start inside the regions;
  clusters     EVAL sfs_subcluster rows of the call provenance (both sides of the weight threshold; w = w_own),
               start inside the regions; one per interval (larger weight), as in tab:res-cascade.
Objects within 50 bp of an SV breakpoint of the V0.5 truth are left out of the classes (2026-10-14).
Class of an object, first match (the SVs are left out: their breakpoints have their own analysis, bp_points.py;
an object near an SV breakpoint takes the class of what else it lies on; until 2026-10-14 an SV class came first,
output kept in classes_with_sv.tsv):
  SNV          within 10 bp of a clonal somatic SNV (MNVs included)
  INDEL_HP     within 10 bp of a homopolymer indel (adds or removes copies of a base that also flanks it)
  INDEL_other  within 10 bp of any other indel
  no_smvar     start outside the smvar v0.3 regions (no small-variant truth there)
  other_trf    within 50 bp of a TRF/VNTR interval (Severus annotation)
  other        the rest
Variants: the clonal smvar v0.3 variants inside both the smvar and the SV regions, one per allele of a multiallelic
record, in minimal representation; covered by an SFS or by a cluster of weight >= 5 that
overlaps [POS - 10, END + 10].

Output: classes.tsv (index, object, class, n, share, median weight), clusters.<index>.tsv.gz (interval, w, class),
        smvar_cov.tsv.gz (per variant), smvar_cov_summary.tsv
Usage: python sfs_classes.py
"""
import bisect, gzip, os, re, subprocess, sys
from collections import Counter, defaultdict
import pysam

PAD_SV, PAD_SM, PAD_TRF, WMIN = 50, 10, 50, 5
A = os.environ["HG008"]
R = f"{A}/reads_ec_align"
O = os.environ["HG008"] + "/reads_ec_align/sfs_signal"; os.makedirs(O, exist_ok=True)
SVREG = f"{A}/ground_truths/v0.5/GRCh38_HG008-T-V0.5_somatic-stvar-clonal_and_subclonal.draftbenchmark.bed"
GT = f"{A}/ground_truths/v0.5/GRCh38_HG008-T-V0.5_somatic-stvar_ALL.draftbenchmark.vcf"
SM = f"{A}/ground_truths/smvar_v0.3/HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz"
SMREG = f"{A}/ground_truths/smvar_v0.3/HG008-T_somatic_smvar_benchmark_v0.3_all.bed"
TRF = f"{A}/refs/sev_extras/human_GRCh38_no_alt_analysis_set.trf.bed"
REF = f"{R}/ref/GRCh38_GIABv3.chr1-22X.fa"
IDX = [("utg", f"{R}/germline_cascade/rows/R2_refutg_nofilters"),
       ("full", None),
       ("v63", f"{R}/germline_cascade/rows/R2w_refv63_nofilters")]
CLASSES = ["SNV", "INDEL_HP", "INDEL_other", "no_smvar", "other_trf", "other"]


def paths(name, d):
    if name == "full": return f"{R}/bp_coverage/sfs_clean.bed.gz", f"{R}/bp_coverage/nogl_w5/provenance.tsv"
    return f"{d}/sfs_clean.bed.gz", f"{d}/provenance.tsv"


def merge(d):
    out = {}
    for c, L in d.items():
        L.sort(); m = []
        for s, e in L:
            if m and s <= m[-1][1]: m[-1][1] = max(m[-1][1], e)
            else: m.append([s, e])
        out[c] = ([a for a, _ in m], [b for _, b in m])
    return out


def hit(t, c, s, e):
    """[s, e) overlaps a merged interval of t."""
    if c not in t: return False
    st, en = t[c]; i = bisect.bisect_right(st, e - 1) - 1
    return i >= 0 and en[i] > s


def inside(t, c, p):
    if c not in t: return False
    st, en = t[c]; i = bisect.bisect_right(st, p) - 1
    return i >= 0 and en[i] > p


def bed(p):
    d = defaultdict(list)
    for l in open(p):
        f = l.split("\t"); d[f[0]].append((int(f[1]), int(f[2])))
    return d


def intersect(a, b):
    out = {}
    for c in set(a) & set(b):
        x, y = a[c], b[c]; i = j = 0; m = []
        while i < len(x[0]) and j < len(y[0]):
            s, e = max(x[0][i], y[0][j]), min(x[1][i], y[1][j])
            if s < e: m.append((s, e))
            if x[1][i] < y[1][j]: i += 1
            else: j += 1
        if m: out[c] = ([s for s, _ in m], [e for _, e in m])
    return out


def minimal(p, ref, alt):
    """Minimal representation of one allele: common suffix trimmed, then common prefix down to one anchor base."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]: ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0] and ref[1] == alt[1]: ref, alt, p = ref[1:], alt[1:], p + 1
    return p, ref, alt


def vtype(ref, alt, fa, c, p):
    if len(ref) == len(alt): return "SNV"
    d = alt[1:] if len(alt) > len(ref) else ref[1:]
    if len(set(d)) == 1:
        b = d[0].upper(); after = fa.fetch(c, p + len(ref), p + len(ref) + 1).upper()
        if ref[0].upper() == b or after == b: return "INDEL_HP"
    return "INDEL_other"


def classify(tabs, c, s, e):
    for k in CLASSES[:-1]:
        if k == "no_smvar":
            if not inside(tabs[k], c, s): return k
        elif hit(tabs[k], c, s, e): return k
    return "other"


def median(v):
    v = sorted(v); return v[len(v) // 2] if v else ""


if __name__ == "__main__":
    smreg = merge(bed(SMREG)); reg = merge(bed(SVREG))
    both = intersect(smreg, reg)
    reg_bp = sum(e - s for c in reg for s, e in zip(*reg[c]))
    print(f"[sfs_classes] regions: {reg_bp / 1e9:.3f} Gb", file=sys.stderr, flush=True)
    fa = pysam.FastaFile(REF)
    variants = []; sm = defaultdict(lambda: defaultdict(list))
    for l in gzip.open(SM, "rt"):
        if l.startswith("#"): continue
        f = l.split("\t"); c, p0 = f[0], int(f[1]) - 1
        if not inside(both, c, p0) or not inside(both, c, p0 + len(f[3]) - 1): continue
        # the v0.3 VCF joins some sites into multiallelic records (19, e.g. ACCCC>A,ACC): one variant per allele
        for alt in f[4].split(","):
            p, ref, alt = minimal(p0, f[3], alt)
            t = vtype(ref, alt, fa, c, p); variants.append((c, p, ref, alt, t))
            sm[t][c].append((p - PAD_SM, p + len(ref) + PAD_SM))
    tabs = {"SNV": merge(sm["SNV"]), "INDEL_HP": merge(sm["INDEL_HP"]),
            "INDEL_other": merge(sm["INDEL_other"])}
    trf = defaultdict(list)
    for c, L in bed(TRF).items(): trf[c] = [(s - PAD_TRF, e + PAD_TRF) for s, e in L]
    tabs["other_trf"] = merge(trf); tabs["no_smvar"] = smreg
    # objects within PAD_SV of an SV breakpoint of the V0.5 truth (all 244 records: POS, and END of DEL/DUP/INV) are left
    # out of the class counts (they belong to the breakpoint analysis); they still count for the small-variant coverage
    svb = defaultdict(list)
    for l in open(GT):
        if l.startswith("#"): continue
        f = l.split("\t"); p = int(f[1]) - 1; svb[f[0]].append((p - PAD_SV, p + PAD_SV + 1))
        m = re.search(r"(?:^|;)END=(\d+)", f[7]); t = re.search(r"SVTYPE=(\w+)", f[7]).group(1)
        if m and t in ("DEL", "DUP", "INV"):
            e = int(m.group(1)) - 1; svb[f[0]].append((e - PAD_SV, e + PAD_SV + 1))
    svt = merge(svb); nsv = {}

    rows = []; vcov = {}
    for name, d in IDX:
        sfs, prov = paths(name, d)
        print(f"[sfs_classes] {name}", file=sys.stderr, flush=True)
        # placed SFS
        cnt = Counter(); sfs_iv = defaultdict(list)
        zc = subprocess.Popen(["zcat", sfs], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for l in zc.stdout:
            f = l.split("\t", 4); c, s, e = f[0], int(f[1]), int(f[2])
            sfs_iv[c].append((s, max(e, s + 1)))
            if inside(reg, c, s):
                if hit(svt, c, s, max(e, s + 1)): cnt["near_sv"] += 1
                else: cnt[classify(tabs, c, s, max(e, s + 1))] += 1
        zc.wait()
        nsv[(name, "SFS")] = cnt.pop("near_sv", 0)
        n = sum(cnt.values())
        for k in CLASSES: rows.append((name, "SFS", k, cnt[k], cnt[k] / n, ""))
        rows.append((name, "SFS", "total", n, 1.0, ""))
        # clusters
        # one cluster per interval, as in tab:res-cascade: the clusters of the two haplotypes of a locus can share
        # the interval; the cluster takes the larger weight
        cw = defaultdict(list); cl_iv = defaultdict(list); iv = {}
        for l in open(prov):
            f = l.rstrip("\n").split("\t")
            if f[0] != "EVAL" or f[1] != "sfs_subcluster": continue
            m = re.match(r"S:([^:]+):(\d+)-(\d+)", f[2]); key = (m.group(1), int(m.group(2)), int(m.group(3)) + 1)
            iv[key] = max(iv.get(key, 0), int(f[3]))
        with gzip.open(f"{O}/clusters.{name}.tsv.gz", "wt") as fo:
            fo.write("chrom\tstart\tend\tw\tclass\n")
            for (c, s, e), w in sorted(iv.items()):
                if w >= WMIN: cl_iv[c].append((s, e))
                if not inside(reg, c, s): continue
                if hit(svt, c, s, e): nsv[(name, "clusters")] = nsv.get((name, "clusters"), 0) + (w >= WMIN); continue
                k = classify(tabs, c, s, e); cw[k].append(w)
                fo.write(f"{c}\t{s}\t{e}\t{w}\t{k}\n")
        for thr, lab in ((WMIN, f"clusters w>={WMIN}"), (0, "clusters all")):
            n = sum(1 for k in CLASSES for w in cw[k] if w >= thr)
            for k in CLASSES:
                ws = [w for w in cw[k] if w >= thr]
                rows.append((name, lab, k, len(ws), len(ws) / n if n else 0, median(ws)))
            rows.append((name, lab, "total", n, 1.0, median([w for k in CLASSES for w in cw[k] if w >= thr])))
        st, ct = merge(sfs_iv), merge(cl_iv)
        for c, p, ref, alt, t in variants:
            s, e = p - PAD_SM, p + len(ref) + PAD_SM
            vcov.setdefault((c, p, ref, alt, t), {})[name] = (int(hit(st, c, s, e)), int(hit(ct, c, s, e)))

    with open(f"{O}/classes.tsv", "w") as fo:
        fo.write(f"# regions: V0.5 clonal_and_subclonal, {reg_bp / 1e9:.3f} Gb\n")
        fo.write("# left out, within 50 bp of an SV breakpoint: " + "; ".join(f"{k[0]} {k[1]} {v}" for k, v in sorted(nsv.items())) + "\n")
        fo.write("index\tobject\tclass\tn\tshare\tmedian_w\n")
        for r in rows: fo.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]:.4f}\t{r[5]}\n")
    names = [n for n, _ in IDX]
    with gzip.open(f"{O}/smvar_cov.tsv.gz", "wt") as fo:
        fo.write("chrom\tpos\tref\talt\ttype\t" + "\t".join(f"{n}_sfs\t{n}_cluster" for n in names) + "\n")
        for (c, p, ref, alt, t), v in vcov.items():
            fo.write(f"{c}\t{p + 1}\t{ref}\t{alt}\t{t}\t" + "\t".join(f"{v[n][0]}\t{v[n][1]}" for n in names) + "\n")
    with open(f"{O}/smvar_cov_summary.tsv", "w") as fo:
        fo.write("type\tn\t" + "\t".join(f"{n}_sfs\t{n}_cluster" for n in names) + "\n")
        for t in ("SNV", "INDEL_HP", "INDEL_other"):
            sel = [v for k, v in vcov.items() if k[4] == t]
            fo.write(f"{t}\t{len(sel)}\t" + "\t".join(
                f"{sum(v[n][0] for v in sel) / len(sel):.4f}\t{sum(v[n][1] for v in sel) / len(sel):.4f}"
                for n in names) + "\n")
    print("done", file=sys.stderr)
