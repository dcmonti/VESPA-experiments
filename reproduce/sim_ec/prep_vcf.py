#!/usr/bin/env python3
"""Prepares a VCF for truvari: (1) PASS only if requested (never on the GT),
(2) the normalizations of bench_events_multicaller.prepare (SVLEN from SVINSLEN for the
<INS> of nanomonsv, END=POS+|SVLEN| for the INV of Severus v1.0), (3) region filter
PER EVENT: a record is kept only if ALL its anchors (POS, END for
DEL/DUP/INV, the mate locus in the ALT for BNDs) fall in the BED. Thus the two
halves of a BND pair stay or leave together, as in bench_events_multicaller
(--regions-bed), instead of truvari's --includebed, which looks at one record at a time.
usage: prep_vcf.py IN.vcf OUT.vcf REGIONS.bed|- [--pass] [--minlen N]   ("-" = no region filter)
--minlen: discards the non-BND records with |SVLEN| (or END-POS) < N. Needed on the GT for truvari:
filter_mate_event.py counts the GT events on the file, not on what truvari kept."""
import sys, os, re
from bisect import bisect_right
from collections import defaultdict
sys.path.insert(0, os.environ["LIB"])
import bench_events_multicaller as bm
src, dst, bed = sys.argv[1:4]; passonly = "--pass" in sys.argv
minlen = int(sys.argv[sys.argv.index("--minlen") + 1]) if "--minlen" in sys.argv else 0
reg = defaultdict(list)
for l in (open(bed) if bed != "-" else []):
    if l.strip() and not l.startswith(("#", "track")):
        c = l.split("\t"); reg[c[0]].append((int(c[1]), int(c[2])))
for k in reg: reg[k].sort()
st = {k: [x[0] for x in v] for k, v in reg.items()}
def inside(c, p):
    v = reg.get(c)
    if not v: return False
    j = bisect_right(st[c], p - 1) - 1
    return j >= 0 and v[j][1] > p - 1
BND = re.compile(r"[\[\]]([^\[\]:]+):(\d+)[\[\]]")
tmp, *_ = bm.prepare(src, passonly)
kept = drop = 0
with open(dst, "w") as out:
    for l in open(tmp):
        if l.startswith("#"): out.write(l); continue
        c = l.rstrip("\n").split("\t")
        anch = [(c[0], int(c[1]))]
        m = BND.search(c[4])
        if m: anch.append((m.group(1), int(m.group(2))))
        else:
            t = re.search(r"(?:^|;)SVTYPE=([A-Za-z]+)", c[7])
            e = re.search(r"(?:^|;)END=(\d+)", c[7])
            if t and t.group(1).upper() not in ("INS", "BND") and e:
                anch.append((c[0], int(e.group(1))))
        if minlen and not m:
            sl = re.search(r"(?:^|;)SVLEN=(-?\d+)", c[7])
            size = abs(int(sl.group(1))) if sl else (anch[1][1] - anch[0][1] if len(anch) > 1 else abs(len(c[4]) - len(c[3])))
            if size < minlen: drop += 1; continue
        if bed == "-" or all(inside(a, p) for a, p in anch): out.write(l); kept += 1
        else: drop += 1
os.unlink(tmp)
# INFO keys used in the records but not declared in the header (e.g. the SAVANA VCFs of
# Qin & Li): bcftools tolerates them, pysam writing (truvari) does not. They are declared
# (END Integer, keys without value Flag, the others String/Number=.); no value changes.
L = open(dst).read().split("\n")
declared = {re.search(r"ID=([^,>]+)", l).group(1) for l in L if l.startswith("##INFO=")}
used, valued = set(), set()
for l in L:
    if l and not l.startswith("#"):
        for kv in l.split("\t")[7].split(";"):
            if kv and kv != ".":
                used.add(kv.split("=")[0])
                if "=" in kv: valued.add(kv.split("=")[0])
used.add("END")  # pysam writes END on symbolic records even if the file does not use it
missing = sorted(used - declared)
if missing:
    def decl(k):
        if k == "END": return f'##INFO=<ID={k},Number=1,Type=Integer,Description="declared by prep_vcf.py">'
        if k not in valued: return f'##INFO=<ID={k},Number=0,Type=Flag,Description="declared by prep_vcf.py">'
        return f'##INFO=<ID={k},Number=.,Type=String,Description="declared by prep_vcf.py">'
    i = next(i for i, l in enumerate(L) if l.startswith("#CHROM"))
    L[i:i] = [decl(k) for k in missing]
    open(dst, "w").write("\n".join(L))
    print(f"{os.path.basename(src)}: declared INFO {missing}", file=sys.stderr)
# Samples with the same name (Sniffles2 multi-sample of Qin & Li: "grch38 grch38"):
# bcftools rejects the header. They are renamed with a suffix; the records do not change.
L = open(dst).read().split("\n")
i = next(i for i, l in enumerate(L) if l.startswith("#CHROM"))
cols = L[i].split("\t"); seen = {}
for k in range(9, len(cols)):
    n = cols[k]; seen[n] = seen.get(n, 0) + 1
    if seen[n] > 1: cols[k] = f"{n}_{seen[n]}"
if "\t".join(cols) != L[i]:
    L[i] = "\t".join(cols); open(dst, "w").write("\n".join(L))
    print(f"{os.path.basename(src)}: renamed duplicate sample names", file=sys.stderr)
import subprocess
subprocess.run(["bcftools", "sort", "-o", dst + ".s", dst], check=True, stderr=subprocess.DEVNULL)
os.replace(dst + ".s", dst)
print(f"{os.path.basename(src)}: kept {kept} dropped {drop}", file=sys.stderr)
