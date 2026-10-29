#!/usr/bin/env python3
"""Converts the minisv tabular format (the "*_filterasm.vcf" files of Qin & Li, without a VCF
header) into a VCF: columns ctg1 pos1 ori ctg2 pos2 . count strand INFO. DEL/DUP/INV/INS
become symbolic records (END = pos2, SVLEN from the INFO or from pos2-pos1); junctions
between different chromosomes (and intra-chromosomal BNDs) become BNDs with the bracket
derived from the orientation: '>>' t[p[, '<<' ]p]t, '><' t]p], '<>' [p[t (convention of the
gsv/minisv formats; to be verified, matters only for truvari, which checks the orientation).
7-column lines are indels from the CIGAR: ctg start end . count strand INFO.
usage: msv2vcf.py IN OUT"""
import sys, re
src, dst = sys.argv[1:3]
ALT = {">>": "N[{c}:{p}[", "<<": "]{c}:{p}]N", "><": "N]{c}:{p}]", "<>": "[{c}:{p}[N"}
hdr = ['##fileformat=VCFv4.2', '##source=minisv (converted by msv2vcf.py)',
       '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
       '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">',
       '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
       '##INFO=<ID=MSV_COUNT,Number=1,Type=Integer,Description="minisv read count (col 7)">',
       '##FILTER=<ID=PASS,Description="All filters passed">']
contigs, recs = [], []
for n, l in enumerate(open(src)):
    c = l.rstrip("\n").split("\t")
    if len(c) == 7:   # indel from the CIGAR: ctg start end . count strand INFO
        c1, p1, ori, c2, p2, cnt, info = c[0], int(c[1]), ">>", c[0], int(c[2]), c[4], c[6]
    else:             # junction: ctg1 pos1 ori ctg2 pos2 . count strand INFO
        c1, p1, ori, c2, p2, cnt, info = c[0], int(c[1]), c[2], c[3], int(c[4]), c[6], c[8]
    t = re.search(r"SVTYPE=(\w+)", info).group(1)
    sl = re.search(r"SVLEN=(-?\d+)", info)
    for x in (c1, c2):
        if x not in contigs: contigs.append(x)
    if c1 == c2 and t in ("DEL", "DUP", "INV", "INS"):
        L = abs(int(sl.group(1))) if sl else abs(p2 - p1)
        end = p1 if t == "INS" else p2
        recs.append((c1, p1, f"msv{n}", "N", f"<{t}>", f"SVTYPE={t};SVLEN={-L if t == 'DEL' else L};END={end};MSV_COUNT={cnt}"))
    else:
        recs.append((c1, p1, f"msv{n}", "N", ALT[ori].format(c=c2, p=p2), f"SVTYPE=BND;MSV_COUNT={cnt}"))
with open(dst, "w") as o:
    o.write("\n".join(hdr + [f"##contig=<ID={x}>" for x in contigs]) + "\n")
    o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for r in recs:
        o.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t.\tPASS\t{r[5]}\n")
print(f"{src}: {len(recs)} records", file=sys.stderr)
