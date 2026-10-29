#!/usr/bin/env python3
"""Known somatic variants on hap1 and tumor haplotypes hap1_T1 (clonal) and hap1_T2 (clonal +
subclonal), with the truth in reference coordinates (hg38 region, sequence CTG).

Positions are chosen only in the blocks where hap1 is aligned 1:1 to the reference
(ref/hap1_vs_ref.paf, long cs): >= 1 kb from every germline indel and >= 50 bp from every germline
SNP, outside the N. Thus the local sequence of hap1 and of the reference is identical and the
position is carried to the reference with a single offset, without ambiguity.
Minimum distances: 500 bp between short variants, 5 kb between an SV and any other variant.
Types (2026-09-29: SNV and short indels only if they change the sequence compressed at cap 2 and if the +-25 bp window with
the alternative allele, compressed, does not already occur in hap1, hap2 or the reference; INS > 500 bp copies only):
  SNV          different base;
  INDEL_HP     homopolymer of >= 4 bases lengthened or shortened by 1 base (the real case; at cap 2
               the run stays >= 2 before and after, so the change is invisible);
  INDEL_other  insertion or deletion of 2-20 bp outside homopolymers (flanks differ from the
               inserted or removed bases, so the representation is unique);
  DEL 50 bp-5 kb, INS 50 bp-2 kb (half random sequence, half copy of a segment of the region),
  INV and tandem DUP 0.5-10 kb (log-uniform lengths).
The variants of each type are assigned in turn to the clonal and subclonal clones.
USAGE: python3 02_somatic.py D CTG SEED N_SNV N_IND N_HP N_DEL N_INS N_INV N_DUP
Output: ref/hap1_T1.fa, ref/hap1_T2.fa, truth/truth_small.vcf, truth/truth_sv.vcf, truth/truth.tsv
"""
import sys, random, re, math, os
D, CTG, SEED = sys.argv[1], sys.argv[2], int(sys.argv[3])
N_SNV, N_IND, N_HP, N_DEL, N_INS, N_INV, N_DUP = map(int, sys.argv[4:11])
rng = random.Random(SEED)
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def read_fa(p):
    s = []
    for l in open(p):
        if not l.startswith(">"): s.append(l.strip())
    return "".join(s).upper()

ref = read_fa(f"{D}/ref/region.fa"); hap = read_fa(f"{D}/ref/hap1.fa")
hap2 = read_fa(f"{D}/ref/hap2.fa") if os.path.exists(f"{D}/ref/hap2.fa") else ""

# --- hap1 -> reference map from the blocks of the long cs ----------------------------------------
# 1:1 blocks (q, t, len); germline differences in hap1 coordinates
blocks, gsnp, gind = [], [], []
for l in open(f"{D}/ref/hap1_vs_ref.paf"):
    f = l.rstrip("\n").split("\t")
    if f[4] != "+" or "tp:A:P" not in f: continue
    cs = next(x[5:] for x in f[12:] if x.startswith("cs:Z:"))
    q, t = int(f[2]), int(f[7])
    for op, val in re.findall(r"([=*+\-~])([A-Za-z]+(?:[a-z])?)", cs):
        if op == "=":
            blocks.append((q, t, len(val))); q += len(val); t += len(val)
        elif op == "*":
            gsnp.append(q); q += 1; t += 1
        elif op == "+":               # bases present in hap1 (query), absent from the reference
            gind.append((q, q + len(val))); q += len(val)
        elif op == "-":               # reference bases absent from hap1
            gind.append((q, q + 1)); t += len(val)
blocks.sort(); gsnp.sort(); gind.sort()
bstart = [b[0] for b in blocks]
import bisect
def to_ref(q):
    """position in hap1 -> position in the reference, if q is in a 1:1 block, otherwise None"""
    i = bisect.bisect_right(bstart, q) - 1
    if i < 0: return None
    bq, bt, bl = blocks[i]
    return bt + (q - bq) if q < bq + bl else None

# "safe" positions: in a block, away from germline indels (1 kb) and SNPs (50 bp), no N
safe = bytearray(len(hap))
for bq, bt, bl in blocks:
    safe[bq:bq + bl] = b"\x01" * bl
def clear(a, b):
    a = max(0, a); b = min(len(hap), b)
    if a < b: safe[a:b] = b"\x00" * (b - a)
for a, b in gind: clear(a - 1000, b + 1000)
for p in gsnp: clear(p - 50, p + 51)
for m in re.finditer(r"N+", hap): clear(m.start() - 1000, m.end() + 1000)
def all_safe(a, b): return b <= len(hap) and a >= 0 and all(safe[a:b])
print(f"# hap1 {len(hap)} bp, ref {len(ref)} bp, 1:1 blocks {len(blocks)}, germline SNPs {len(gsnp)}, "
      f"germline indels {len(gind)}, safe bases {sum(safe)}", file=sys.stderr)

# --- choice of the positions with the minimum distances ------------------------------------------
taken = []   # (a, b, is_sv)
def free(a, b, sv):
    pad = 5000 if sv else 500
    for x, y, s in taken:
        p = 5000 if (s or sv) else 500
        if a < y + p and b > x - p: return False
    return True
def loglen(lo, hi): return int(round(math.exp(rng.uniform(math.log(lo), math.log(hi)))))
def rand_seq(n): return "".join(rng.choice("ACGT") for _ in range(n))
def hp_ok(seq, a, b):
    """no homopolymer >= 3 crossing the boundaries [a, b)"""
    w = seq[max(0, a - 3):b + 3]
    return not re.search(r"A{3,}|C{3,}|G{3,}|T{3,}", w)

variants = []   # dict: id, type, clone, q (hap1, 0-based), apply info, truth info
def pick(kind, n, make, sv, tries=200000):
    k = 0
    for _ in range(tries):
        if k == n: break
        v = make()
        if v is None: continue
        a, b = v["span"]
        # short variants: the whole interval +-50 bp safe; SV: only the ends (+-50 bp), the interior
        # may contain germline differences and the end is carried to the reference with the map
        ok = (all_safe(a - 50, a + 50) and all_safe(b - 50, b + 50)) if sv else all_safe(a - 50, b + 50)
        if not ok or not free(a, b, sv): continue
        taken.append((a, b, sv))
        v["type"] = kind; v["clone"] = "truncal" if k % 2 == 0 else "subclonal"
        variants.append(v); k += 1
    if k < n: sys.exit(f"{kind}: found only {k} positions of {n}, relax the constraints or change region")

L = len(hap)
# SVs first (bulkier), then the short variants
def mk_del():
    n = loglen(50, 5000); q = rng.randrange(1, L - n - 1)
    return {"q": q, "len": n, "span": (q, q + n)}
def mk_ins():
    n = loglen(50, 2000); q = rng.randrange(1, L - 1)
    # above 500 bp copies only (like mobile elements): a long novel sequence produces, after
    # assembly of the overlapping SFS, a single SFS as long as the insertion, removed by the cut at
    # 500 bp; in the real data (HG008 V0.5) the INS > 500 bp are all ME (2026-09-29)
    if n <= 500 and rng.random() < 0.5: s = rand_seq(n); src = "random"
    else:
        a = rng.randrange(0, L - n); s = hap[a:a + n]; src = f"copy:{a}"
        if "N" in s: return None
    return {"q": q, "len": n, "seq": s, "src": src, "span": (q, q + 1)}
def mk_inv():
    n = loglen(500, 10000); q = rng.randrange(1, L - n - 1)
    return {"q": q, "len": n, "span": (q, q + n)}
def mk_dup():
    n = loglen(500, 10000); q = rng.randrange(1, L - n - 1)
    return {"q": q, "len": n, "span": (q, q + n)}
def hpc(x, cap=2): return re.sub(r"A+|C+|G+|T+", lambda m: m.group()[0] * min(len(m.group()), cap), x)
def hpc_visible(q, ref_len, alt):
    """does the variant change the sequence compressed at cap 2 in a +-20 bp window? (2026-09-29: the
    variants invisible after HPC would measure HPC and not EC, like the homopolymer indels)"""
    a, b = max(0, q - 20), min(L, q + ref_len + 20)
    return hpc(hap[a:b]) != hpc(hap[a:q] + alt + hap[q + ref_len:b])
# genome of the normal compressed at cap 2 like the index, both strands: a variant whose window with
# the alternative allele already occurs in it cannot produce SFS (2026-09-29: "variants identical to the alt")
NORMAL_HPC = "|".join(hpc(x) + "|" + hpc(x.translate(COMP)[::-1]) for x in (hap, hap2, ref) if x)
def alt_novel(q, ref_len, alt, w=25):
    a, b = max(0, q - w), min(L, q + ref_len + w)
    return hpc(hap[a:q] + alt + hap[q + ref_len:b]) not in NORMAL_HPC
def mk_snv():
    q = rng.randrange(1, L - 1); b = hap[q]
    if b not in "ACGT": return None
    alt = rng.choice([x for x in "ACGT" if x != b])
    if not hpc_visible(q, 1, alt) or not alt_novel(q, 1, alt): return None
    return {"q": q, "alt": alt, "span": (q, q + 1)}
def mk_ind():
    n = rng.randint(2, 20); q = rng.randrange(20, L - n - 20)
    if rng.random() < 0.5:           # deletion of hap[q:q+n]
        if not hp_ok(hap, q, q + n) or hap[q - 1] == hap[q] or hap[q + n - 1] == hap[q + n]: return None
        if not hpc_visible(q, n, "") or not alt_novel(q, n, ""): return None
        return {"q": q, "len": n, "kind": "del", "span": (q, q + n)}
    s = rand_seq(n)                  # insertion between q-1 and q
    if re.search(r"A{3,}|C{3,}|G{3,}|T{3,}", hap[q - 3:q] + s + hap[q:q + 3]) or s[0] == hap[q - 1] or s[-1] == hap[q]: return None
    if not hpc_visible(q, 0, s) or not alt_novel(q, 0, s): return None
    return {"q": q, "len": n, "kind": "ins", "seq": s, "span": (q, q + 1)}
hp_runs = [(m.start(), m.end()) for m in re.finditer(r"A{4,}|C{4,}|G{4,}|T{4,}", hap)]
def mk_hp():
    a, b = rng.choice(hp_runs); n = 1; base = hap[a]
    if rng.random() < 0.5:
        return {"q": a, "len": n, "kind": "del", "base": base, "run": (a, b), "span": (a, b)}
    return {"q": a, "len": n, "kind": "ins", "seq": base * n, "base": base, "run": (a, b), "span": (a, b)}

pick("DEL", N_DEL, mk_del, True); pick("INS", N_INS, mk_ins, True)
pick("INV", N_INV, mk_inv, True); pick("DUP", N_DUP, mk_dup, True)
pick("SNV", N_SNV, mk_snv, False); pick("INDEL_other", N_IND, mk_ind, False); pick("INDEL_HP", N_HP, mk_hp, False)

# --- truth in reference coordinates ---------------------------------------------------------------
for i, v in enumerate(sorted(variants, key=lambda x: x["q"])):
    v["id"] = f"SIM_{v['type']}_{i}"
    q = v["q"]; t = to_ref(q); assert t is not None
    if v["type"] == "SNV":
        assert ref[t] == hap[q]; v["vcf"] = (t + 1, ref[t], v["alt"], {})
    elif v["type"] in ("INDEL_other", "INDEL_HP"):
        if v["type"] == "INDEL_HP":    # left-normalized representation: anchor before the run
            a, b = v["run"]; ta = to_ref(a); assert ref[ta:ta + (b - a)] == hap[a:b]
            if v["kind"] == "del": v["vcf"] = (ta, ref[ta - 1:ta + v["len"]], ref[ta - 1], {})
            else: v["vcf"] = (ta, ref[ta - 1], ref[ta - 1] + v["seq"], {})
        elif v["kind"] == "del":
            assert ref[t:t + v["len"]] == hap[q:q + v["len"]]
            v["vcf"] = (t, ref[t - 1:t + v["len"]], ref[t - 1], {})
        else: v["vcf"] = (t, ref[t - 1], ref[t - 1] + v["seq"], {})
    else:
        n = v["len"]
        if v["type"] == "INS": v["vcf"] = (t, ref[t - 1], "<INS>", {"END": t, "SVLEN": n}); continue
        te = to_ref(q + n); assert te is not None; m = te - t     # length on the reference
        if v["type"] == "DEL": v["vcf"] = (t, ref[t - 1], "<DEL>", {"END": te, "SVLEN": -m})
        elif v["type"] == "INV": v["vcf"] = (t + 1, ref[t], "<INV>", {"END": te, "SVLEN": m})
        else: v["vcf"] = (t + 1, ref[t], "<DUP>", {"END": te, "SVLEN": m})

# --- tumor haplotypes: edits applied from right to left -----------------------------------------
def apply(vs):
    s = list(hap)
    for v in sorted(vs, key=lambda x: -x["q"]):
        q, ty = v["q"], v["type"]
        if ty == "SNV": s[q] = v["alt"]
        elif ty in ("INDEL_other", "INDEL_HP"):
            if v["kind"] == "del": del s[q:q + v["len"]]
            else: s[q:q] = list(v["seq"])
        elif ty == "DEL": del s[q:q + v["len"]]
        elif ty == "INS": s[q:q] = list(v["seq"])
        elif ty == "INV": s[q:q + v["len"]] = list("".join(s[q:q + v["len"]]).translate(COMP)[::-1])
        elif ty == "DUP": s[q + v["len"]:q + v["len"]] = s[q:q + v["len"]]
    return "".join(s)
def write_fa(p, name, seq):
    with open(p, "w") as o:
        o.write(f">{name}\n")
        for i in range(0, len(seq), 80): o.write(seq[i:i + 80] + "\n")
tr = [v for v in variants if v["clone"] == "truncal"]
write_fa(f"{D}/ref/hap1_T1.fa", "hap1_T1", apply(tr))
write_fa(f"{D}/ref/hap1_T2.fa", "hap1_T2", apply(variants))

os.makedirs(f"{D}/truth", exist_ok=True)
hdr = ["##fileformat=VCFv4.2", f"##contig=<ID={CTG},length={len(ref)}>",
       '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">',
       '##INFO=<ID=END,Number=1,Type=Integer,Description="end">',
       '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="length">',
       '##INFO=<ID=CLONE,Number=1,Type=String,Description="truncal (hap1_T1 and hap1_T2) or subclonal (hap1_T2 only)">',
       '##INFO=<ID=VTYPE,Number=1,Type=String,Description="SNV, INDEL_HP, INDEL_other, DEL, INS, INV, DUP">',
       '##ALT=<ID=DEL,Description="deletion">', '##ALT=<ID=INS,Description="insertion">',
       '##ALT=<ID=INV,Description="inversion">', '##ALT=<ID=DUP,Description="tandem duplication">',
       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
small = sorted([v for v in variants if v["type"] in ("SNV", "INDEL_other", "INDEL_HP")], key=lambda v: v["vcf"][0])
svs = sorted([v for v in variants if v["type"] in ("DEL", "INS", "INV", "DUP")], key=lambda v: v["vcf"][0])
for name, vs in (("truth_small.vcf", small), ("truth_sv.vcf", svs)):
    with open(f"{D}/truth/{name}", "w") as o:
        o.write("\n".join(hdr) + "\n")
        for v in vs:
            pos, r, a, extra = v["vcf"]
            info = [f"SVTYPE={v['type']}" if v["type"] in ("DEL", "INS", "INV", "DUP") else None,
                    *(f"{k}={x}" for k, x in extra.items()), f"CLONE={v['clone']}", f"VTYPE={v['type']}"]
            o.write(f"{CTG}\t{pos}\t{v['id']}\t{r}\t{a}\t.\tPASS\t{';'.join(x for x in info if x)}\n")
with open(f"{D}/truth/truth.tsv", "w") as o:
    o.write("id\ttype\tclone\thap1_pos0\tref_pos\tlen\tdetail\n")
    for v in sorted(variants, key=lambda x: x["q"]):
        o.write(f"{v['id']}\t{v['type']}\t{v['clone']}\t{v['q']}\t{v['vcf'][0]}\t{v.get('len', 1)}\t{v.get('src', v.get('kind', '.'))}\n")
from collections import Counter
print("# variants:", dict(Counter((v["type"], v["clone"]) for v in variants)), file=sys.stderr)
