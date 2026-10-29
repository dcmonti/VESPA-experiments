#!/usr/bin/env python3
"""Prepare a VCF (truth or calls) for Minda and truvari.

1. --pass-only: keep records with FILTER PASS or '.' (never used on the truth).
2. Representation fixes, no change of content:
   - INS with SVINSLEN but no SVLEN (nanomonsv): SVLEN = SVINSLEN;
   - non-INS, non-BND records with END = POS and a non-zero SVLEN (Severus 1.0 inversions): END = POS + |SVLEN|.
3. --min-len N: drop non-BND records shorter than N bp (|SVLEN|, else END - POS).
4. --bed: keep a record only if all its anchors fall in the BED: POS, END for DEL/DUP/INV, and the mate
   locus in the ALT of a BND. The two records of a BND pair are therefore kept or dropped together.
5. FILTER is set to PASS on every record kept (Minda reads only PASS records).
6. INFO keys used but not declared in the header are declared (pysam, used by truvari, rejects them),
   duplicate sample names get a suffix (bcftools rejects them), and the output is sorted.

Usage: prep_vcf.py IN.vcf OUT.vcf [--bed REGIONS.bed] [--pass-only] [--min-len N]
"""
import argparse
import os
import re
import subprocess
import sys
from bisect import bisect_right
from collections import defaultdict

BND_ALT = re.compile(r"[\[\]]([^\[\]:]+):(\d+)[\[\]]")
SVTYPE = re.compile(r"(?:^|;)SVTYPE=([A-Za-z]+)")
END = re.compile(r"(?:^|;)END=(\d+)")
SVLEN = re.compile(r"(?:^|;)SVLEN=(-?\d+)")
SVINSLEN = re.compile(r"(?:^|;)SVINSLEN=(-?\d+)")


def read_bed(path):
    regions = defaultdict(list)
    for line in open(path):
        if line.strip() and not line.startswith(("#", "track")):
            c = line.split("\t")
            regions[c[0]].append((int(c[1]), int(c[2])))
    for chrom in regions:
        regions[chrom].sort()
    starts = {chrom: [s for s, _ in v] for chrom, v in regions.items()}
    return regions, starts


def inside(regions, starts, chrom, pos):
    """1-based position pos in a 0-based half-open BED interval."""
    v = regions.get(chrom)
    if not v:
        return False
    j = bisect_right(starts[chrom], pos - 1) - 1
    return j >= 0 and v[j][1] > pos - 1


def fix_record(c):
    """Representation fixes of step 2, on the split record c."""
    info = c[7]
    if "SVTYPE=INS" in info and not SVLEN.search(info):
        m = SVINSLEN.search(info)
        if m:
            info = f"SVLEN={m.group(1)};{info}"
    t = SVTYPE.search(info)
    if t and t.group(1).upper() not in ("INS", "BND") and "[" not in c[4] and "]" not in c[4]:
        e, sl = END.search(info), SVLEN.search(info)
        if e and sl and int(e.group(1)) == int(c[1]) and abs(int(sl.group(1))) > 0:
            info = info[:e.start(1)] + str(int(c[1]) + abs(int(sl.group(1)))) + info[e.end(1):]
    c[7] = info
    return c


def anchors(c):
    """Positions that must fall in the BED, and whether the record is a BND."""
    a = [(c[0], int(c[1]))]
    m = BND_ALT.search(c[4])
    if m:
        a.append((m.group(1), int(m.group(2))))
        return a, True
    t, e = SVTYPE.search(c[7]), END.search(c[7])
    if t and t.group(1).upper() not in ("INS", "BND") and e:
        a.append((c[0], int(e.group(1))))
    return a, False


def size(c, a):
    sl = SVLEN.search(c[7])
    if sl:
        return abs(int(sl.group(1)))
    if len(a) > 1:
        return a[1][1] - a[0][1]
    return abs(len(c[4]) - len(c[3]))


def fix_header(lines):
    """Declare INFO keys used but not declared; rename duplicate sample names."""
    declared = {re.search(r"ID=([^,>]+)", l).group(1) for l in lines if l.startswith("##INFO=")}
    used, valued = {"END"}, set()  # pysam writes END on symbolic records even when the file does not use it
    for l in lines:
        if l and not l.startswith("#"):
            for kv in l.split("\t")[7].split(";"):
                if kv and kv != ".":
                    used.add(kv.split("=")[0])
                    if "=" in kv:
                        valued.add(kv.split("=")[0])
    new = []
    for k in sorted(used - declared):
        if k == "END":
            new.append(f'##INFO=<ID={k},Number=1,Type=Integer,Description="declared by prep_vcf.py">')
        elif k not in valued:
            new.append(f'##INFO=<ID={k},Number=0,Type=Flag,Description="declared by prep_vcf.py">')
        else:
            new.append(f'##INFO=<ID={k},Number=.,Type=String,Description="declared by prep_vcf.py">')
    i = next(i for i, l in enumerate(lines) if l.startswith("#CHROM"))
    cols, seen = lines[i].split("\t"), {}
    for k in range(9, len(cols)):
        seen[cols[k]] = seen.get(cols[k], 0) + 1
        if seen[cols[k]] > 1:
            cols[k] = f"{cols[k]}_{seen[cols[k]]}"
    lines[i] = "\t".join(cols)
    lines[i:i] = new
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vcf_in")
    ap.add_argument("vcf_out")
    ap.add_argument("--bed", help="keep only records with all anchors in these regions")
    ap.add_argument("--pass-only", action="store_true", help="keep only FILTER PASS or '.'")
    ap.add_argument("--min-len", type=int, default=0, help="drop non-BND records shorter than this")
    args = ap.parse_args()

    regions, starts = read_bed(args.bed) if args.bed else (None, None)
    lines, kept, dropped = [], 0, 0
    for line in open(args.vcf_in, errors="replace"):
        line = line.rstrip("\n")
        if line.startswith("#"):
            lines.append(line)
            continue
        c = line.split("\t")
        if len(c) < 8 or (args.pass_only and c[6] not in ("PASS", ".")):
            continue
        c = fix_record(c)
        a, is_bnd = anchors(c)
        if (args.min_len and not is_bnd and size(c, a) < args.min_len) or \
                (regions is not None and not all(inside(regions, starts, ch, p) for ch, p in a)):
            dropped += 1
            continue
        c[6] = "PASS"
        lines.append("\t".join(c))
        kept += 1

    with open(args.vcf_out + ".unsorted", "w") as out:
        out.write("\n".join(fix_header(lines)) + "\n")
    subprocess.run(["bcftools", "sort", "-o", args.vcf_out, args.vcf_out + ".unsorted"],
                   check=True, stderr=subprocess.DEVNULL)
    os.remove(args.vcf_out + ".unsorted")
    print(f"{os.path.basename(args.vcf_in)}: kept {kept}, dropped {dropped}", file=sys.stderr)


if __name__ == "__main__":
    main()
