#!/usr/bin/env python3
"""Split the svviz2 input of COLO829 (colo829_gt_svviz2.vcf) into N chunk VCFs for the events that have no report yet
in reports/; a BND pair (same EVENT, which svviz2 uses as the variant name) is always in one chunk. 2026-10-27."""
import glob, os, re
N = 8
C = os.environ["COLO829"] + "/reads_ec_align_s50/svviz2_vaf"
done = {os.path.basename(p).split(".")[0] for p in glob.glob(f"{C}/reports/*.report.tsv")}
header, groups = [], {}
for l in open(f"{C}/colo829_gt_svviz2.vcf"):
    if l.startswith("#"): header.append(l); continue
    c = l.split("\t")
    name = re.search(r"EVENT=([^;\t]*)", c[7]).group(1) if "SVTYPE=BND" in c[7] else c[2]
    groups.setdefault(name, []).append(l)
todo = sorted(g for g in groups if g not in done)
for k in range(N):
    with open(f"{C}/chunks/chunk_{k}.vcf", "w") as f:
        f.writelines(header)
        rows = [l for g in todo[k::N] for l in groups[g]]
        f.writelines(sorted(rows, key=lambda l: (l.split("\t")[0], int(l.split("\t")[1]), l.split("\t")[2])))
print(f"variants {len(groups)}, done {len(done & set(groups))}, to run {len(todo)}")
