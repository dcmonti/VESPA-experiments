#!/usr/bin/env python3
"""Project SFS (read coordinates) onto the reference via the tumor EC BAM.

Usage:
    sfs_to_genome.py <sfs> <bam> <out.bed.gz> [threads] [--primary-only]

Output (TSV, gzipped): chrom, start, end (0-based half-open), read, qstart,
sfs_len, exact, aln.

  exact = 1  the SFS overlaps an aligned block, so the position is real
  exact = 0  the SFS lies entirely in a clipped tail of the PRIMARY alignment
             and no alignment places it: the position is snapped to the
             alignment boundary and is only indicative
  aln   = P  primary alignment, S supplementary

**Supplementary alignments are used.**  A SFS spanning an SV junction sits in
the clipped tail of the primary and is properly placed only by the
supplementary record covering that part of the read -- exactly the case that
matters for BND and large DUP.  Such a SFS therefore gets one row per
alignment that covers it (dedup on read+qstart if you need one row per SFS).

Coordinates: the SFS are forward-strand (search was run on FASTA), so they are
flipped for reverse alignments; hard clips count towards the read-level query
offset (minimap2 writes H in supplementary records).
"""
import gzip
import sys
from multiprocessing import Pool

import pysam


def load_sfs(path):
    sfs = {}
    last = ""
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            r, s, l = c[0], int(c[1]), int(c[2])
            if r == "*":
                r = last
            else:
                last = r
            sfs.setdefault(r, []).append((s, l))
    return sfs


SFS = {}
BAM = ""
PRIMARY_ONLY = False


def aligned_blocks(aln):
    """[(q_start, q_end, ref_start)] for the aligned (M/=/X) blocks, with q in
    whole-read coordinates in the orientation of this record (H and S count)."""
    blocks = []
    q = 0
    ref = aln.reference_start
    for op, ln in aln.cigartuples:
        if op in (4, 5):            # S, H: query only, not aligned
            q += ln
        elif op == 1:               # I: query only
            q += ln
        elif op in (2, 3):          # D, N: reference only
            ref += ln
        elif op in (0, 7, 8):       # M, =, X
            blocks.append((q, q + ln, ref))
            q += ln
            ref += ln
    return blocks


def place(blocks, a, b):
    """Reference span of the query interval [a, b] (inclusive), or None when it
    falls entirely outside the aligned blocks."""
    lo = hi = None
    for qs, qe, rs in blocks:
        if b < qs or a >= qe:
            continue
        s = max(a, qs)
        e = min(b, qe - 1)
        r1, r2 = rs + (s - qs), rs + (e - qs)
        if lo is None or r1 < lo:
            lo = r1
        if hi is None or r2 > hi:
            hi = r2
    return None if lo is None else (lo, hi)


def work(chrom):
    out = []
    bam = pysam.AlignmentFile(BAM, "rb", threads=2)
    for aln in bam.fetch(chrom):
        if aln.is_secondary or aln.is_unmapped:
            continue
        if PRIMARY_ONLY and aln.is_supplementary:
            continue
        hits = SFS.get(aln.query_name)
        if not hits:
            continue
        L = aln.infer_read_length()
        if not L:
            continue
        blocks = aligned_blocks(aln)
        if not blocks:
            continue
        tag = "S" if aln.is_supplementary else "P"
        for s, l in hits:
            e = s + l - 1
            a, b = (L - 1 - e, L - 1 - s) if aln.is_reverse else (s, e)
            if a < 0 or b >= L:
                continue
            span = place(blocks, a, b)
            if span is not None:
                out.append((chrom, span[0], span[1] + 1, aln.query_name,
                            s, l, 1, tag))
            elif not aln.is_supplementary:
                # entirely inside a clipped tail: snap to the nearest boundary,
                # flagged as not exact.  A supplementary elsewhere may still
                # place it properly, hence the dedup step downstream.
                ref = blocks[0][2] if b < blocks[0][0] else blocks[-1][2] + (
                    blocks[-1][1] - blocks[-1][0]) - 1
                out.append((chrom, ref, ref + 1, aln.query_name, s, l, 0, "P"))
    bam.close()
    return out


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    sfs_path, bam_path, out_path = args[0], args[1], args[2]
    nthreads = int(args[3]) if len(args) > 3 else 16

    global SFS, BAM, PRIMARY_ONLY
    PRIMARY_ONLY = "--primary-only" in sys.argv

    bam = pysam.AlignmentFile(bam_path, "rb")
    chroms = sorted(zip(bam.references, bam.lengths), key=lambda x: -x[1])
    bam.close()

    print("loading SFS..", file=sys.stderr, flush=True)
    SFS = load_sfs(sfs_path)
    BAM = bam_path
    print(f"loaded SFS for {len(SFS)} reads", file=sys.stderr, flush=True)

    n = exact = supp = 0
    with gzip.open(out_path, "wt") as fo, Pool(nthreads) as pool:
        for rows in pool.imap_unordered(work, [c for c, _ in chroms],
                                        chunksize=1):
            for r in rows:
                fo.write("\t".join(map(str, r)) + "\n")
                exact += r[6]
                supp += r[7] == "S"
            n += len(rows)
            print(f"[{n}] done chunk of {len(rows)}", file=sys.stderr,
                  flush=True)
    print(f"placed {n} rows ({exact} exact, {supp} from supplementary)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
