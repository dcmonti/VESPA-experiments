#!/usr/bin/env python3
"""
Converts a VCF in which intra-chromosomal SVs are encoded as BREAKEND PAIRS
into a VCF with symbolic records (<DEL>/<DUP>/<INV>/<INS>), which is the representation
that truvari can compare.

WHY IT IS NEEDED. `gt/truthset_COLO829_hg38lifted.vcf` is the Valle-Inclan et al.
2022 truth set (68 somatic SVs of COLO829) lifted over to GRCh38. The file encodes *everything* in
breakend notation, deletions included:

    1  117973563  truthset_2_1  T  T[1:117973599[  SVLEN=36;SVTYPE=DEL;MATEID=truthset_2_2
    1  117973599  truthset_2_2  T  ]1:117973563]T  SVLEN=36;SVTYPE=DEL;MATEID=truthset_2_1

132 of the 134 records are like this. Truvari compares DEL/INS/DUP/INV by position and
size on a symbolic or sequence-resolved ALT: a DEL written as a BND pair
cannot match a `<DEL>` of VESPA, so the bench gave TP=1 (only the 2 symbolic
records in the file). It is not a liftover problem: it is a
representation problem. This script solves it.

WHAT IT DOES
  - groups records by event (ID `<base>_1` / `<base>_2`, cross-checked with MATEID)
  - SVTYPE in DEL/DUP/INV/INS with both mates on the SAME chromosome -> a single
    symbolic record: POS=min(p1,p2), END=max(p1,p2), ALT=<SVTYPE>,
    SVLEN negative for DEL and positive for the others (convention of VESPA and of the
    Zenodo v2 benchmark)
  - SVTYPE=BND -> the two lines stay unchanged (they are inter-chromosomal translocations,
    and VESPA emits them as BND: the representation is already comparable)
  - single records already symbolic -> passed through, with END added if missing
  - adds the `chr` prefix to CHROM, to the header contigs and to the chromosomes inside the
    breakend notation of the ALTs (the original file uses `1`, our BAMs/VCFs use `chr1`)

SVLEN AND LIFTOVER. In 3 events the SVLEN in the file (computed on GRCh37) does not match
the distance between the breakends lifted to GRCh38 (truthset_22: 24534 vs 32532; truthset_24:
67479 vs 67480; truthset_62: 49 vs 43). The GRCh38 positions are the authoritative
information in the lifted file, so SVLEN is RECOMPUTED from END-POS; the original
is kept in INFO/SVLEN_ORIG and the event is marked with INFO/LEN_MISMATCH, so that
those 3 events can be excluded from the evaluation. No other record is touched.

PRESERVED METADATA. INFO/SUPP_SEQ (supporting sequencing technologies),
INFO/SUPP_VAL (validation: PCR / CAPTURE / BIONANO / NOT_VALIDATED), INFO/GENE,
INFO/CLUSTER pass through unchanged: SUPP_VAL in particular allows weighting the
confidence of individual events at bench time (52 of the 68 events are PCR-validated,
2 are NOT_VALIDATED).

USAGE
  python3 lib/bnd_pairs_to_symbolic.py IN.vcf -o OUT.vcf [--ref hg38.fa[.gz]]
  # --ref: re-reads the REF base from the reference at the new POS (recommended: the REF in
  #        the file comes from a 2021 liftover and may no longer match). Without
  #        --ref the REF of the source record is kept.
  # --no-add-chr: do not add the chr prefix (default: added if absent)
"""

import argparse
import re
import sys
from collections import OrderedDict, defaultdict

SYMBOLIC = ("DEL", "DUP", "INV", "INS")
# SVLEN negative for deletions, positive for the rest (VESPA / Zenodo v2 convention)
NEGATIVE_SVLEN = ("DEL",)

BND_ALT = re.compile(r"^(?P<pre>[A-Za-z.]*)(?P<br1>[\[\]])(?P<chrom>[^:\[\]]+):(?P<pos>\d+)(?P<br2>[\[\]])(?P<post>[A-Za-z.]*)$")
EVENT_ID = re.compile(r"^(?P<base>.+?)_(?P<idx>\d+)$")


def parse_info(s):
    d = OrderedDict()
    if s == ".":
        return d
    for field in s.split(";"):
        if not field:
            continue
        k, _, v = field.partition("=")
        d[k] = v if _ else None
    return d


def fmt_info(d):
    out = []
    for k, v in d.items():
        out.append(k if v is None else "%s=%s" % (k, v))
    return ";".join(out) if out else "."


class Rec(object):
    __slots__ = ("chrom", "pos", "vid", "ref", "alt", "qual", "filt", "info", "rest")

    def __init__(self, line):
        f = line.rstrip("\n").split("\t")
        self.chrom, self.pos, self.vid, self.ref, self.alt, self.qual, self.filt = (
            f[0], int(f[1]), f[2], f[3], f[4], f[5], f[6])
        self.info = parse_info(f[7])
        self.rest = f[8:]

    def line(self):
        f = [self.chrom, str(self.pos), self.vid, self.ref, self.alt,
             self.qual, self.filt, fmt_info(self.info)] + list(self.rest)
        return "\t".join(f)


def add_chr(name):
    return name if name.startswith("chr") else "chr" + name


def fix_alt_chrom(alt):
    """Adds chr to the chromosome inside an ALT in breakend notation."""
    m = BND_ALT.match(alt)
    if not m:
        return alt
    g = m.groupdict()
    return "%s%s%s:%s%s%s" % (g["pre"], g["br1"], add_chr(g["chrom"]), g["pos"],
                              g["br2"], g["post"])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vcf")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--ref", help="reference FASTA: re-reads the REF base at the new POS")
    ap.add_argument("--no-add-chr", action="store_true",
                    help="do not add the chr prefix")
    args = ap.parse_args()

    fasta = None
    if args.ref:
        try:
            import pysam
        except ImportError:
            sys.exit("--ref requires pysam (conda activate sv_somatic)")
        fasta = pysam.FastaFile(args.ref)

    meta, chrom_line, recs = [], None, []
    with open(args.vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                meta.append(line.rstrip("\n"))
            elif line.startswith("#CHROM"):
                chrom_line = line.rstrip("\n")
            elif line.strip():
                recs.append(Rec(line))
    if chrom_line is None:
        sys.exit("VCF without #CHROM line")

    want_chr = not args.no_add_chr and not all(r.chrom.startswith("chr") for r in recs)

    # ---- group by event ----
    groups = defaultdict(list)
    order = []
    for r in recs:
        m = EVENT_ID.match(r.vid)
        base = m.group("base") if m else r.vid
        if base not in groups:
            order.append(base)
        groups[base].append(r)

    stats = defaultdict(int)
    warnings = []
    out_recs = []

    for base in order:
        grp = groups[base]
        svtype = grp[0].info.get("SVTYPE", ".")

        # consistency of the pairing declared in MATEID
        if len(grp) == 2:
            ids = {grp[0].vid, grp[1].vid}
            mates = {grp[0].info.get("MATEID"), grp[1].info.get("MATEID")}
            if None not in mates and mates != ids:
                warnings.append("%s: inconsistent MATEID (%s vs %s)"
                                % (base, sorted(mates), sorted(ids)))

        # --- translocations: left as a breakend pair ---
        if svtype == "BND":
            for r in grp:
                if want_chr:
                    r.chrom = add_chr(r.chrom)
                    r.alt = fix_alt_chrom(r.alt)
                    if "MATEID" in r.info:
                        pass  # the IDs do not contain chromosomes
                r.info["EVENT"] = base
                out_recs.append(r)
            stats["bnd_kept"] += 1
            continue

        # --- single record, already symbolic ---
        if len(grp) == 1:
            r = grp[0]
            if not r.alt.startswith("<"):
                warnings.append("%s: single non-symbolic record (ALT=%s), left"
                                " unchanged" % (base, r.alt))
                stats["passthrough"] += 1
            else:
                # INFO rewritten in the same canonical order as the collapsed records,
                # so the two groups can be queried with the same expressions
                svlen = abs(int(float(r.info.get("SVLEN", 0))))
                info = OrderedDict()
                info["SVTYPE"] = svtype
                info["SVLEN"] = str(-svlen if svtype in NEGATIVE_SVLEN else svlen)
                info["END"] = str(r.pos if svtype == "INS" else r.pos + svlen)
                info["CHR2"] = add_chr(r.chrom) if want_chr else r.chrom
                info["EVENT"] = base
                for k in ("SUPP_SEQ", "SUPP_VAL", "GENE", "CLUSTER"):
                    if k in r.info:
                        info[k] = r.info[k]
                r.info = info
                stats["single_symbolic"] += 1
            if want_chr:
                r.chrom = add_chr(r.chrom)
            r.vid = base
            r.info.setdefault("EVENT", base)
            out_recs.append(r)
            continue

        if len(grp) != 2:
            warnings.append("%s: %d records, expected 1 or 2 -- left unchanged"
                            % (base, len(grp)))
            for r in grp:
                if want_chr:
                    r.chrom = add_chr(r.chrom)
                    r.alt = fix_alt_chrom(r.alt)
                out_recs.append(r)
            stats["passthrough"] += len(grp)
            continue

        a, b = sorted(grp, key=lambda r: r.pos)
        if a.chrom != b.chrom:
            warnings.append("%s: SVTYPE=%s but mates on different chromosomes (%s/%s) --"
                            " left as breakends" % (base, svtype, a.chrom, b.chrom))
            for r in grp:
                if want_chr:
                    r.chrom = add_chr(r.chrom)
                    r.alt = fix_alt_chrom(r.alt)
                r.info["EVENT"] = base
                out_recs.append(r)
            stats["passthrough"] += 2
            continue

        if svtype not in SYMBOLIC:
            warnings.append("%s: SVTYPE=%s not convertible -- left unchanged"
                            % (base, svtype))
            for r in grp:
                if want_chr:
                    r.chrom = add_chr(r.chrom)
                    r.alt = fix_alt_chrom(r.alt)
                out_recs.append(r)
            stats["passthrough"] += 2
            continue

        # --- collapse of the pair into a symbolic record ---
        chrom = add_chr(a.chrom) if want_chr else a.chrom
        span = b.pos - a.pos
        orig = a.info.get("SVLEN")
        orig_abs = abs(int(float(orig))) if orig not in (None, ".", "NA") else None

        info = OrderedDict()
        info["SVTYPE"] = svtype
        if svtype == "INS":
            # for an INS the length is that of the inserted sequence, not the span
            length = orig_abs if orig_abs is not None else span
            info["SVLEN"] = str(length)
            info["END"] = str(a.pos)
        else:
            info["SVLEN"] = str(-span if svtype in NEGATIVE_SVLEN else span)
            info["END"] = str(b.pos)
            if orig_abs is not None and orig_abs != span:
                info["SVLEN_ORIG"] = str(orig_abs)
                info["LEN_MISMATCH"] = None
                stats["len_mismatch"] += 1
                warnings.append("%s: file SVLEN=%d but END-POS=%d -- using END-POS,"
                                " original in SVLEN_ORIG, event marked LEN_MISMATCH"
                                % (base, orig_abs, span))
        info["CHR2"] = chrom
        info["EVENT"] = base
        for k in ("SUPP_SEQ", "SUPP_VAL", "GENE", "CLUSTER"):
            if k in a.info:
                info[k] = a.info[k]

        ref = a.ref
        if fasta is not None:
            try:
                ref = fasta.fetch(chrom, a.pos - 1, a.pos).upper() or a.ref
            except (KeyError, ValueError) as exc:
                warnings.append("%s: REF not readable from %s (%s), keeping '%s'"
                                % (base, chrom, exc, a.ref))
        if len(ref) != 1:
            ref = ref[:1] or "N"

        new = Rec("\t".join([chrom, str(a.pos), base, ref, "<%s>" % svtype,
                             a.qual, a.filt, "."] + list(a.rest)))
        new.info = info
        out_recs.append(new)
        stats["collapsed_" + svtype] += 1

    # ---- header ----
    drop = ("##contig", "##ALT", "##INFO", "##fileformat")
    kept = [m for m in meta if not m.startswith(drop)]
    contigs = [m for m in meta if m.startswith("##contig")]
    if want_chr:
        contigs = [re.sub(r"(<ID=)([^,>]+)", lambda m: m.group(1) + add_chr(m.group(2)),
                          c) for c in contigs]

    header = ["##fileformat=VCFv4.2"]
    header += kept
    header.append("##source=bnd_pairs_to_symbolic.py (breakend pairs collapsed into"
                  " symbolic records; SVLEN recomputed from END-POS)")
    header += contigs
    for t, d in (("DEL", "Deletion"), ("DUP", "Duplication"), ("INV", "Inversion"),
                 ("INS", "Insertion")):
        header.append('##ALT=<ID=%s,Description="%s">' % (t, d))
    header += [
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV; negative for deletions. Recomputed as END-POS for collapsed breakend pairs">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the SV">',
        '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome of the second breakend">',
        '##INFO=<ID=EVENT,Number=1,Type=String,Description="Truth-set event ID shared by the two breakends of an SV">',
        '##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of the mate breakend">',
        '##INFO=<ID=SVLEN_ORIG,Number=1,Type=Integer,Description="SVLEN as given in the input VCF, kept when it disagrees with END-POS (GRCh37-era length preserved through liftover)">',
        '##INFO=<ID=LEN_MISMATCH,Number=0,Type=Flag,Description="Input SVLEN disagrees with the distance between the lifted breakends">',
        '##INFO=<ID=SUPP_SEQ,Number=.,Type=String,Description="Supporting sequencing technologies">',
        '##INFO=<ID=SUPP_VAL,Number=.,Type=String,Description="Supporting validation technologies (PCR, CAPTURE, BIONANO, NOT_VALIDATED)">',
        '##INFO=<ID=GENE,Number=.,Type=String,Description="Gene(s) overlapping the breakpoint">',
        '##INFO=<ID=CLUSTER,Number=.,Type=String,Description="IDs of other clustered breakpoints">',
    ]

    # sorting by the order of the header contigs
    corder = {}
    for c in contigs:
        m = re.search(r"<ID=([^,>]+)", c)
        if m:
            corder[m.group(1)] = len(corder)
    out_recs.sort(key=lambda r: (corder.get(r.chrom, len(corder)), r.pos))

    with open(args.out, "w") as fh:
        fh.write("\n".join(header) + "\n")
        fh.write(chrom_line + "\n")
        for r in out_recs:
            fh.write(r.line() + "\n")

    # ---- report ----
    sys.stderr.write("Input events: %d (%d records)\n" % (len(order), len(recs)))
    for k in sorted(stats):
        sys.stderr.write("  %-22s %d\n" % (k, stats[k]))
    sys.stderr.write("Output records: %d -> %s\n" % (len(out_recs), args.out))
    if warnings:
        sys.stderr.write("\nWarnings (%d):\n" % len(warnings))
        for w in warnings:
            sys.stderr.write("  - %s\n" % w)


if __name__ == "__main__":
    main()
