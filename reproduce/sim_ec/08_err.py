#!/usr/bin/env python3
"""Residual error of the tumor reads of one cell relative to the sequence they were simulated from
(2026-09-30; possible only on simulated data). Each read (cells/<S>/<EC>/tumor.ec.fa.gz: corrected,
or raw for EC=raw) goes to its own sequence, from the name <sample>.<component>.<n>: hap1, hap2, hap1_T1,
hap1_T2 (ref/<component>.fa, with the somatic variants already applied). Thus neither the germline differences
from hg38 nor the somatic variants enter: what remains is the sequencing error left (or introduced) by
EC; a somatic allele erased by EC counts as an error.
minimap2 -c -x map-hifi --eqx, primary alignments only (tp:A:P); from the CIGAR: X, bases of I and D.
Writes cells/<S>/<EC>/eval/err.tsv. USAGE: python3 08_err.py D S EC [THREADS]"""
import sys, os, re, gzip, shutil, subprocess
D, S, EC = sys.argv[1:4]; TH = sys.argv[4] if len(sys.argv) > 4 else "16"
W = f"{D}/cells/{S}/{EC}"; E = f"{W}/eval"; TMP = f"{W}/err_tmp"
os.makedirs(E, exist_ok=True); os.makedirs(TMP, exist_ok=True)
MM2 = os.environ["MINIMAP2"]
COMPS = ("hap1", "hap2", "hap1_T1", "hap1_T2")
comp = lambda r: r.split(".")[-2]

# --- reads split by component ---------------------------------------------------------------------
out = {c: open(f"{TMP}/{c}.fa", "w") for c in COMPS}; nreads = {c: 0 for c in COMPS}; nbases = {c: 0 for c in COMPS}
name = None
for l in gzip.open(f"{W}/tumor.ec.fa.gz", "rt"):
    if l[0] == ">": name = l[1:].split()[0]; c = comp(name); out[c].write(l); nreads[c] += 1
    else: out[c].write(l); nbases[c] += len(l.rstrip("\n"))
for f in out.values(): f.close()

# --- alignment to their own sequence ----------------------------------------------------------------
CG = re.compile(r"(\d+)([=XID])")
tot = {c: [0, 0, 0, 0, 0] for c in COMPS}   # aligned reads, aligned bases (query), X, I, D
for c in COMPS:
    if not nreads[c]: continue
    p = subprocess.Popen([MM2, "-c", "-x", "map-hifi", "--eqx", "--secondary=no", "-t", TH, f"{D}/ref/{c}.fa", f"{TMP}/{c}.fa"],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    for l in p.stdout:
        f = l.rstrip("\n").split("\t")
        if "tp:A:P" not in f[12:]: continue
        cg = next(x[5:] for x in f[12:] if x.startswith("cg:Z:"))
        t = tot[c]; t[0] += 1; t[1] += int(f[3]) - int(f[2])
        for n, op in CG.findall(cg):
            if op == "X": t[2] += int(n)
            elif op == "I": t[3] += int(n)
            elif op == "D": t[4] += int(n)
    if p.wait(): sys.exit(f"minimap2 failed on {c}")
shutil.rmtree(TMP)

with open(f"{E}/err.tsv", "w") as o:
    o.write("sample\tec\tcomponent\treads\tbases\treads_aligned\tbases_aligned\tmismatch\tins_bp\tdel_bp\terr_per_kb\n")
    rows = [(c, *tot[c]) for c in COMPS if nreads[c]]
    allt = [sum(r[i] for r in rows) for i in range(1, 6)]
    for c, ra, ba, x, i, d in rows + [("all", *allt)]:
        nr = sum(nreads.values()) if c == "all" else nreads[c]; nb = sum(nbases.values()) if c == "all" else nbases[c]
        o.write(f"{S}\t{EC}\t{c}\t{nr}\t{nb}\t{ra}\t{ba}\t{x}\t{i}\t{d}\t{1000 * (x + i + d) / ba if ba else float('nan'):.4f}\n")
print(f"{S} {EC}: residual error {1000 * sum(allt[2:]) / allt[1]:.3f}/kb over {allt[1]} aligned bases ({E}/err.tsv)")
