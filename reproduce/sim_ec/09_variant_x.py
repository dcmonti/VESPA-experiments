#!/usr/bin/env python3
"""Per-variant test of the error-correction model (final_version.tex, sec:ec) on the simulation.

The model: a somatic variant is carried by X reads, X ~ NB(mu = C rho f / 2, alpha); the correction keeps a base
supported by at least T reads, so a variant is lost when X < T. Here X is known read by read, because the name of
every read gives the sequence it was simulated from (<sample>.<component>.<n>): hap1_T1 and hap1_T2 carry the clonal
variants, hap1_T2 only the subclonal ones.

Per cell (sample x correction setting) and short variant of truth/truth_small.vcf:
  x_origin   reads of the carrying components whose primary alignment spans the variant (SNV: the base; indel: POS-WIN
             to POS+len(REF)+WIN), in the raw BAM;
  x_alt_raw  of those, the reads that show the alternative allele in the raw BAM (a sequencing error at the base
             removes a read): the support the correction sees;
  x_alt_ec   of the x_alt_raw reads, those that still show it in the corrected BAM (same read names);
  n_alt_ec   all reads showing it in the corrected BAM, any origin;
  sfs        covered at the SFS stage with the rule of 07_attrib.py (an SFS of prov/sfs_clean.bed.gz containing the
             variant, from a read showing the alternative allele in the corrected BAM); reproduces eval/attr.tsv;
  sfs_origin the same, restricted to the x_alt_raw reads.
Allele in a read: SNV, the read base at POS equals ALT; indel, an indel of the same type and length starting within
+-WIN bp of POS (07_attrib.py, alt_reads). All alignments of the simulation are primary with MAPQ 60 (checked
2026-10-14 on t60_p0.3: 16 942 primary reads, none below MAPQ 60, none unmapped, same reads before and after the
correction), so no mapping-quality filter is applied.

Output:
  cells/<S>/<EC>/eval/variant_x.tsv           one row per variant;
  table/variant_x_summary.tsv                 per cell, type and clone: n, T, mu of the model (hifiasm coverage C of
                                              the NB log x rho x f / 2), mean and variance of x_origin and x_alt_raw,
                                              alpha by moments of each, fraction with x_alt_raw < T, fraction lost (no
                                              read keeps the allele; no SFS);
  table/variant_x_by_x.tsv                    per cell, type, clone and x_alt_raw: n, lost (reads), lost (SFS).
USAGE: python3 09_variant_x.py [D=$SIM_DIR or sim_ec_v2_2026-09-29] [PROCS=15] [CELLS=all, e.g. t60_p0.3/fnr0.10]
"""
import sys, os, re, gzip, glob, csv, statistics
from collections import defaultdict
from multiprocessing import Pool
import pysam

A = os.environ["HG008"]
D = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] else os.environ.get("SIM_DIR", f"{A}/sim_ec_v2_2026-09-29")
PROCS = int(sys.argv[2]) if len(sys.argv) > 2 else 15
ONLY = set(sys.argv[3].split(",")) if len(sys.argv) > 3 else None
CTG, WIN = "sim", 5
CARRY = {"truncal": {"hap1_T1", "hap1_T2"}, "subclonal": {"hap1_T2"}}
CLONE = {"truncal": "clonal", "subclonal": "subclonal"}
F = {"truncal": 1.0, "subclonal": 0.5}            # fraction of tumor cells carrying the variant
comp = lambda r: r.split(".")[-2]


def truth():
    out = []
    for l in open(f"{D}/truth/truth_small.vcf"):
        if l[0] == "#": continue
        f = l.rstrip("\n").split("\t")
        info = dict(x.split("=") for x in f[7].split(";") if "=" in x)
        out.append((f[2], info["VTYPE"], info["CLONE"], int(f[1]) - 1, f[3], f[4]))
    return out


def has_alt(a, p, ref, alt):
    """True if the alignment a shows the alternative allele of the variant at p (0-based)"""
    snv = len(ref) == len(alt) == 1
    ins = len(alt) > len(ref); L = abs(len(alt) - len(ref))
    rp, qp = a.reference_start, 0
    for op, n in a.cigartuples:
        if snv:
            if op in (0, 7, 8):
                if rp <= p < rp + n: return a.query_sequence[qp + p - rp].upper() == alt.upper()
                rp += n; qp += n
            elif op in (2, 3):
                if rp <= p < rp + n: return False
                rp += n
            elif op in (1, 4): qp += n
            if rp > p: return False
        else:
            if op in (0, 7, 8): rp += n
            elif op == 1:
                if ins and n == L and abs(rp - (p + 1)) <= WIN: return True
            elif op == 2:
                if not ins and n == L and abs(rp - (p + 1)) <= WIN: return True
                rp += n
            elif op == 3: rp += n
            if rp > p + WIN + len(ref): return False
    return False


def window(p, ref, alt):
    return (p, p + 1) if len(ref) == len(alt) == 1 else (p - WIN, p + len(ref) + WIN)


def reads_at(bam, p, ref, alt):
    """primary alignments spanning the variant window, by name"""
    s, e = window(p, ref, alt); out = {}
    for a in bam.fetch(CTG, max(0, s), e):
        if a.is_secondary or a.is_supplementary or a.is_unmapped: continue
        if a.reference_start <= s and a.reference_end >= e: out[a.query_name] = a
    return out


def threshold(W):
    """T used by hifiasm (NB cells) or None; hifiasm coverage C of the NB log"""
    # the NB line is in hifiasm.log (threshold.txt keeps only the last lines of the log)
    t = open(f"{W}/hifiasm.log", errors="ignore").read() if os.path.exists(f"{W}/hifiasm.log") else ""
    m = re.search(r"cov=(\d+) .*=> het_cov=(\d+)", t)
    return (int(m.group(2)), int(m.group(1))) if m else (None, None)


def cell(job):
    S, EC = job; W = f"{D}/cells/{S}/{EC}"
    by_read = defaultdict(list)
    for l in gzip.open(f"{W}/prov/sfs_clean.bed.gz", "rt"):
        f = l.split("\t")
        if f[0] == CTG: by_read[f[3]].append((int(f[1]), int(f[2])))
    T, _ = threshold(W)
    raw, ec = pysam.AlignmentFile(f"{D}/raw/{S}.bam"), pysam.AlignmentFile(f"{W}/tumor.ec.bam")
    rows = []
    for vid, vt, cl, p, ref, alt in truth():
        s, e = window(p, ref, alt)
        r_raw, r_ec = reads_at(raw, p, ref, alt), reads_at(ec, p, ref, alt)
        orig = {n for n in r_raw if comp(n) in CARRY[cl]}
        alt_raw = {n for n in orig if has_alt(r_raw[n], p, ref, alt)}
        # 07_attrib.py looks for the allele in every primary alignment overlapping the window, not only in those that
        # span it: same here for the corrected BAM, so that sfs reproduces eval/attr.tsv
        alt_ec_all = set()
        for a in ec.fetch(CTG, max(0, s), e):
            if a.is_secondary or a.is_supplementary or a.is_unmapped: continue
            if has_alt(a, p, ref, alt): alt_ec_all.add(a.query_name)
        alt_ec = alt_raw & alt_ec_all
        cov = lambda rs: any(a < e and b > s for r in rs for a, b in by_read.get(r, ()))
        rows.append((vid, vt, CLONE[cl], p + 1, len(orig), len(alt_raw), len(alt_ec), len(alt_ec_all),
                     int(cov(alt_ec_all)), int(cov(alt_ec))))
    os.makedirs(f"{W}/eval", exist_ok=True)
    with open(f"{W}/eval/variant_x.tsv", "w") as o:
        o.write("id\ttype\tclone\tpos\tx_origin\tx_alt_raw\tx_alt_ec\tn_alt_ec\tsfs\tsfs_origin\n")
        for r in rows: o.write("\t".join(map(str, r)) + "\n")
    return S, EC, T, rows


if __name__ == "__main__":
    jobs = sorted((c.split("/")[-2], c.split("/")[-1]) for c in glob.glob(f"{D}/cells/*/*") if os.path.isdir(c))
    if ONLY: jobs = [j for j in jobs if f"{j[0]}/{j[1]}" in ONLY]
    with Pool(min(PROCS, len(jobs))) as pool: res = pool.map(cell, jobs)
    # hifiasm coverage per sample, from the NB cells (classic and raw cells have no NB log)
    C = {}
    for S, EC in jobs:
        _, c = threshold(f"{D}/cells/{S}/{EC}")
        if c: C.setdefault(S, c)
    for S in {j[0] for j in jobs} - set(C):
        for W in glob.glob(f"{D}/cells/{S}/fnr*"):
            _, c = threshold(W)
            if c: C[S] = c; break
    suf = "" if ONLY is None else ".test"
    mom = lambda v: ((statistics.pvariance(v) - statistics.mean(v)) / statistics.mean(v) ** 2) if statistics.mean(v) > 0 else float("nan")
    with open(f"{D}/table/variant_x_summary{suf}.tsv", "w") as o, open(f"{D}/table/variant_x_by_x{suf}.tsv", "w") as ob:
        o.write("sample\tec\ttype\tclone\tn\tT\tC_hifiasm\tmu_model\tmean_x_origin\tvar_x_origin\talpha_x_origin\t"
                "mean_x_alt_raw\tvar_x_alt_raw\talpha_x_alt_raw\tfrac_x_alt_raw_lt_T\tlost_reads\tlost_sfs\n")
        ob.write("sample\tec\ttype\tclone\tx_alt_raw\tn\tlost_reads\tlost_sfs\n")
        for S, EC, T, rows in sorted(res, key=lambda r: (r[0], r[1])):
            rho = float(S.split("_p")[1])
            for vt in ("SNV", "INDEL_other"):
                for cl, f in (("clonal", 1.0), ("subclonal", 0.5)):
                    sel = [r for r in rows if r[1] == vt and r[2] == cl]
                    if not sel: continue
                    n = len(sel); xo = [r[4] for r in sel]; xa = [r[5] for r in sel]
                    mu = C[S] * rho * f / 2 if S in C else float("nan")
                    lt = f"{sum(1 for x in xa if x < T) / n:.4f}" if T else "NA"
                    o.write(f"{S}\t{EC}\t{vt}\t{cl}\t{n}\t{T if T else 'NA'}\t{C.get(S, 'NA')}\t{mu:.2f}\t"
                            f"{statistics.mean(xo):.2f}\t{statistics.pvariance(xo):.2f}\t{mom(xo):.4f}\t"
                            f"{statistics.mean(xa):.2f}\t{statistics.pvariance(xa):.2f}\t{mom(xa):.4f}\t{lt}\t"
                            f"{sum(1 for r in sel if r[6] == 0) / n:.4f}\t{sum(1 for r in sel if r[8] == 0) / n:.4f}\n")
                    byx = defaultdict(lambda: [0, 0, 0])
                    for r in sel:
                        b = byx[r[5]]; b[0] += 1; b[1] += r[6] == 0; b[2] += r[8] == 0
                    for x in sorted(byx):
                        ob.write(f"{S}\t{EC}\t{vt}\t{cl}\t{x}\t" + "\t".join(map(str, byx[x])) + "\n")
    print(f"{len(res)} cells done", file=sys.stderr)
