#!/usr/bin/env python3
"""Final pass of the HG008 svviz2 run: the groups of svviz2_input.vcf (a BND pair or one record) without any report
in reports/ are written to N chunk VCFs, missing/chunk_<k>.vcf, BND pairs kept together. 2026-10-27."""
import glob, os, re
N = 8
O = os.environ["HG008"] + "/reads_ec_align/svviz2_vaf"
done = {os.path.basename(p).split(".")[0] for p in glob.glob(f"{O}/reports/*.report.tsv")}
header, rec = [], {}
for l in open(f"{O}/svviz2_input.vcf"):
    if l.startswith("#"): header.append(l); continue
    rec[l.split("\t")[2]] = l
mate = lambda l: (re.search(r"MATEID=([^;\t]*)", l) or [None, None])[1]
groups, seen = [], set()
for i in sorted(rec):
    if i in seen: continue
    m = mate(rec[i]) if "SVTYPE=BND" in rec[i] else None
    g = [i] + ([m] if m in rec else []); seen.update(g); groups.append(g)
todo = [g for g in groups if not done & set(g)]
os.makedirs(f"{O}/missing", exist_ok=True)
for k in range(N):
    rows = [rec[i] for g in todo[k::N] for i in g]
    with open(f"{O}/missing/chunk_{k}.vcf", "w") as f:
        f.writelines(header); f.writelines(sorted(rows, key=lambda l: (l.split("\t")[0], int(l.split("\t")[1]), l.split("\t")[2])))
print(f"groups {len(groups)}, with report {len(groups) - len(todo)}, to run {len(todo)}:", " ".join("/".join(g) for g in todo))
