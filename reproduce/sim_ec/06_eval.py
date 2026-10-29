#!/usr/bin/env python3
"""Evaluation of one cell (sample x EC setting) against the simulated truth.
Writes to cells/<S>/<EC>/eval/:
  small.tsv   short variants by type (SNV, INDEL_HP, INDEL_other) and clone: n, median of the alternative
              reads, with >= MIN_ALT alternative reads after EC (primary, MAPQ 60), with an SFS within
              PAD bp, with a cluster of any weight, with a cluster of weight >= 5 (provenance dump);
  bp.tsv      SV breakpoints by stage (SFS, cluster of any weight, cluster weight >= 5,
              candidates) and by clone, tolerance 50 bp (bp_point_table_prov.py);
  calls.tsv   calls with the normal at weights WEIGHTS, Minda (--tolerance 500 --min_size 50), TP/FP/FN
              totals and recall per clone (from the missed events, <sample>_<caller>_fn.tsv);
  noise.tsv   SFS, clusters of weight >= 5, mismatches/kb of the corrected reads, threshold used by hifiasm.
Types and alternative reads: functions of het6/hpc_caps/sfs_strat/small_variant_support.py.
USAGE: python3 06_eval.py D S EC "WEIGHTS" [PAD=10]"""
import sys, os, re, gzip, subprocess, importlib.util, csv
from collections import defaultdict
D, S, EC, WEIGHTS = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split()
PAD = int(sys.argv[5]) if len(sys.argv) > 5 else 10
A = os.environ["HG008"]
SRC = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("svs", f"{SRC}/small_variant_support.py")
sys.argv = [sys.argv[0], "1", str(PAD)]
svs = importlib.util.module_from_spec(spec); spec.loader.exec_module(svs)
import pysam
W = f"{D}/cells/{S}/{EC}"; E = f"{W}/eval"; os.makedirs(E, exist_ok=True)
CTG = "sim"; REF = f"{D}/ref/region.fa"
ENV = dict(os.environ, PATH=os.environ["SV_BIN"] + ":" + os.environ["PATH"])
TOK = re.compile(r"S:([^:|>+\t]+):(\d+)-(\d+)")

# --- objects of the cell -------------------------------------------------------------------------
sf, cl_any, cl5 = defaultdict(list), defaultdict(list), defaultdict(list)
for l in gzip.open(f"{W}/prov/sfs_clean.bed.gz", "rt"):
    f = l.split("\t"); sf[f[0]].append((int(f[1]), int(f[2])))
n_cl5 = set()
for l in open(f"{W}/prov/provenance.tsv"):           # m_min 5: clusters that pass the threshold
    f = l.rstrip("\n").split("\t")
    if f[0] != "EVAL" or f[-1] != "1": continue
    for c, a, b in TOK.findall(f[2]): cl5[c].append((int(a), int(b) + 1))
    n_cl5.add(f[2])
for l in open(f"{W}/prov_w1/provenance.tsv"):        # m_min 1: all clusters, of any weight
    f = l.rstrip("\n").split("\t")
    if f[0] != "EVAL": continue
    for c, a, b in TOK.findall(f[2]): cl_any[c].append((int(a), int(b) + 1))
T_sf, T_any, T_5 = svs.tree(sf), svs.tree(cl_any), svs.tree(cl5)

# --- short variants ------------------------------------------------------------------------------
fa = pysam.FastaFile(REF); rows = []
raw = pysam.AlignmentFile(f"{D}/raw/{S}.bam")          # raw reads of the sample: alt before EC
with pysam.AlignmentFile(f"{W}/tumor.ec.bam") as bam:
    for l in open(f"{D}/truth/truth_small.vcf"):
        if l[0] == "#": continue
        f = l.rstrip("\n").split("\t"); p = int(f[1]) - 1; ref, alt = f[3], f[4]
        info = dict(x.split("=") for x in f[7].split(";") if "=" in x)
        t = svs.vtype(ref, alt, fa, f[0], p)
        dp, na = svs.count_alt(bam, f[0], p, ref, alt, t)
        _, nr = svs.count_alt(raw, f[0], p, ref, alt, t)
        s, e = p - PAD, p + len(ref) + PAD
        rows.append((info["VTYPE"], info["CLONE"], nr, na,
                     int(svs.hit(T_sf, f[0], s, e)), int(svs.hit(T_any, f[0], s, e)), int(svs.hit(T_5, f[0], s, e))))
with open(f"{E}/small.tsv", "w") as o:
    o.write("sample\tec\ttype\tclone\tn\traw_median_alt\traw_in_reads\tmedian_alt\tin_reads\twith_SFS\twith_cluster_any\twith_cluster_w5\n")
    for ty in ("SNV", "INDEL_HP", "INDEL_other"):
        for cl in ("truncal", "subclonal"):
            x = [r for r in rows if r[0] == ty and r[1] == cl]
            if not x: continue
            ar = sorted(r[2] for r in x); a = sorted(r[3] for r in x); n = len(x)
            o.write(f"{S}\t{EC}\t{ty}\t{cl}\t{n}\t{ar[n//2]}\t{sum(r[2] >= svs.MIN_ALT for r in x)}\t"
                    f"{a[n//2]}\t{sum(r[3] >= svs.MIN_ALT for r in x)}\t"
                    f"{sum(r[4] for r in x)}\t{sum(r[5] for r in x)}\t{sum(r[6] for r in x)}\n")

# --- SV breakpoints by stage and by clone ----------------------------------------------------------
BPTAB = f"{SRC}/bp_point_table_prov.py"
with open(f"{E}/bp.tsv", "w") as o:
    o.write("sample\tec\tclone\tstage\tn_objects\tpoints\tpoints_cov\trecall\n")
    for cl in ("truncal", "subclonal"):
        g = f"{E}/truth_sv.{cl}.vcf"
        with open(g, "w") as fo:
            for l in open(f"{D}/truth/truth_sv.vcf"):
                if l[0] == "#" or f"CLONE={cl}" in l: fo.write(l)
        out = subprocess.run(["python3", BPTAB, f"{W}/prov", f"{W}/prov/provenance.tsv", S], capture_output=True, text=True,
                             env=dict(ENV, GT=g, FAI=f"{REF}.fai", CHROMS=CTG, TOL="50", REGIONS="-",
                                      PROV_ANY=f"{W}/prov_w1/provenance.tsv"), check=True).stdout
        for l in out.splitlines():
            if l.startswith("#") or l.startswith("label"): continue
            f = l.split("\t")
            o.write(f"{S}\t{EC}\t{cl}\t{f[2]}\t{f[3]}\t{f[5]}\t{f[6]}\t{f[7]}\n")

# --- calls with the normal: Minda per weight, recall per clone -----------------------------------
PREP = f"{SRC}/prep_vcf.py"
MINDA = os.environ["MINDA"]; PYM = os.environ["MINDA_PYTHON"]
M = f"{E}/minda"; os.makedirs(f"{M}/vcf", exist_ok=True)
clone = {l.split("\t")[2]: re.search(r"CLONE=(\w+)", l).group(1) for l in open(f"{D}/truth/truth_sv.vcf") if l[0] != "#"}
subprocess.run(["python", PREP, f"{D}/truth/truth_sv.vcf", f"{M}/vcf/gt.vcf", "-"], env=ENV, check=True, stderr=subprocess.DEVNULL)
with open(f"{M}/callers.tsv", "w") as ct:
    for w in WEIGHTS:
        subprocess.run(["python", PREP, f"{W}/call_w{w}/calls.vcf", f"{M}/vcf/w{w}.vcf", "-", "--pass"], env=ENV, check=True, stderr=subprocess.DEVNULL)
        ct.write(f"{M}/vcf/w{w}.vcf\tw{w}\n")
subprocess.run([PYM, MINDA, "truthset", "--base", f"{M}/vcf/gt.vcf", "--tsv", f"{M}/callers.tsv", "--out_dir", f"{M}/out",
                "--tolerance", "500", "--min_size", "50", "--sample_name", "SIM"], env=ENV, check=True,
               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
tot = defaultdict(int)
for v in clone.values(): tot[v] += 1
res = {r["Caller"]: r for r in csv.DictReader(open(f"{M}/out/results/overall.tsv"), delimiter="\t")}
with open(f"{E}/calls.tsv", "w") as o:
    o.write("sample\tec\tweight\tTP\tFP\tFN\tprecision\trecall\tF1\trecall_truncal\trecall_subclonal\n")
    for w in WEIGHTS:
        r = res[f"w{w}"]; missed = defaultdict(int)
        fn = f"{M}/out/SIM_w{w}_fn.tsv"
        if os.path.exists(fn):
            for l in open(fn):
                if l.startswith("#"): continue
                i = l.split("\t")[2]
                if i in clone: missed[clone[i]] += 1
        rt = 1 - missed["truncal"] / tot["truncal"] if tot["truncal"] else float("nan")
        rs = 1 - missed["subclonal"] / tot["subclonal"] if tot["subclonal"] else float("nan")
        o.write(f"{S}\t{EC}\t{w}\t{r['True Positives']}\t{r['False Positives']}\t{r['False Negatives']}\t"
                f"{float(r['Precision']):.4f}\t{float(r['Recall']):.4f}\t{float(r['F1 Score']):.4f}\t{rt:.4f}\t{rs:.4f}\n")

# --- noise and threshold used ----------------------------------------------------------------------
n_sfs = sum(1 for l in open(f"{W}/tumor.dedup.cap500.sfs") if l[0] != "#")
mm = bm = 0
for l in open(f"{W}/stats.txt"):
    f = l.rstrip("\n").split("\t")
    if f[0] == "SN" and f[1] == "mismatches:": mm = int(f[2])
    if f[0] == "SN" and f[1] == "bases mapped (cigar):": bm = int(f[2])
thr = open(f"{W}/threshold.txt").read().strip().splitlines()
het = next((re.search(r"het_cov=(\d+)", x).group(1) for x in thr if "NB somatic threshold" in x), None)
if het is None: het = next((re.search(r"peak_het: (-?\d+)", x).group(1) for x in thr if "peak_het" in x), ".")
with open(f"{E}/noise.tsv", "w") as o:
    o.write("sample\tec\thet_cov_used\tn_sfs\tn_clusters_w5\tmismatch_per_kb\n")
    o.write(f"{S}\t{EC}\t{het}\t{n_sfs}\t{len(n_cl5)}\t{1000 * mm / bm if bm else float('nan'):.4f}\n")
print(f"{S} {EC}: evaluation written to {E}")
