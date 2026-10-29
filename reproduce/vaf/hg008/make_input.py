#!/usr/bin/env python3
"""svviz2 input from the HG008 V0.5 ground truth (ALL VCF), restricted to the events of the clonal and subclonal
benchmark (bench/clonal_and_subclonal/out/HG008_support.tsv, base == True) and their BND mates; 2026-10-27.
Normalised as the COLO829 input (colo829/reads_ec_align_s50/svviz2_vaf/make_input.sh):
  - INFO EVENT removed (Number=A in this header: pysam returns a tuple and svviz2 uses EVENT as the file name)
  - <DUP> -> <DUP:TANDEM>, the only duplication svviz2 accepts
  - BND: the sequence inside ALT next to the bracket is set to REF (svviz2 stops on any other "complex event");
    the inserted bases at a junction are thus dropped
  - BND without its mate in the file are dropped (svviz2 genotypes mate pairs)
  - sequence-resolved records: END set to POS + len(REF) - 1 (svviz2 rejects e.g. SV_13, an INS with END at the
    end of the duplicated copy)
Writes svviz2_input.vcf and chunks/chunk_<k>.vcf (N chunks, BND pairs kept together)."""
import csv, os, re, sys
N = 8
G = os.environ["HG008"] + "/ground_truths/v0.5/GRCh38_HG008-T-V0.5_somatic-stvar_ALL.draftbenchmark.vcf"
SUP = os.environ["HG008"] + "/reads_ec_align/bench/clonal_and_subclonal/out/HG008_support.tsv"
O = os.environ["HG008"] + "/reads_ec_align/svviz2_vaf"
header, recs = [], {}
for l in open(G):
    if l.startswith("#"): header.append(l); continue
    c = l.rstrip("\n").split("\t"); recs[c[2]] = c
ids = set()
for r in csv.DictReader(open(SUP), delimiter="\t"):
    if r["base"] == "True": ids.update(x for x in (r["ID_x"], r["ID_y"]) if x in recs)
mate = lambda c: (re.search(r"MATEID=([^;]*)", c[7]) or [None, None])[1]
for i in list(ids):
    m = mate(recs[i])
    if m and m in recs: ids.add(m)
out, dropped, groups, seen = {}, [], [], set()
for i in sorted(ids, key=lambda x: (recs[x][0], int(recs[x][1]), x)):
    c = list(recs[i]); info = c[7]
    c[7] = re.sub(r";?EVENT=[^;]*", "", info).lstrip(";")
    t = re.search(r"SVTYPE=([^;]*)", info).group(1)
    if t == "DUP" and c[4].startswith("<DUP"): c[4] = "<DUP:TANDEM>"
    if t == "BND":
        m = mate(recs[i])
        if not m or m not in ids: dropped.append((i, "no mate")); continue
        a = c[4]
        if re.match(r"^[\[\]][^\[\]]*[\[\]]", a): c[4] = re.match(r"^[\[\]][^\[\]]*[\[\]]", a).group(0) + c[3][-1]
        else: c[4] = c[3][0] + re.search(r"[\[\]][^\[\]]*[\[\]]$", a).group(0)
        c[3] = c[3][0] if not a.startswith(("[", "]")) else c[3][-1]
    if re.fullmatch(r"[ACGTN]+", c[3]) and re.fullmatch(r"[ACGTN]+", c[4]):
        c[7] = re.sub(r"(^|;)END=\d+", rf"\g<1>END={int(c[1]) + len(c[3]) - 1}", c[7])
    out[i] = c
for i in out:
    if i in seen: continue
    m = mate(recs[i]) if "SVTYPE=BND" in out[i][7] else None
    g = [i] + ([m] if m and m in out else [])
    seen.update(g); groups.append(g)
def write(path, members):
    rows = sorted((out[i] for i in members), key=lambda c: (c[0], int(c[1]), c[2]))
    with open(path, "w") as f:
        f.writelines(header); f.writelines("\t".join(c) + "\n" for c in rows)
write(f"{O}/svviz2_input.vcf", list(out))
for k in range(N):
    write(f"{O}/chunks/chunk_{k}.vcf", [i for g in groups[k::N] for i in g])
print(f"records {len(out)}, svviz2 variants {len(groups)}, dropped {dropped}")
