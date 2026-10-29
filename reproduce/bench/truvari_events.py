#!/usr/bin/env python3
"""Event-level TP, FP and FN from a truvari bench output directory.

truvari counts records; a translocation or any event written as two BND mate records counts twice.
Here the two mates of a pair count as one event, on both sides:
  - truth pairs: MATEID/MATE_ID, and BND records whose ALTs name each other's locus (within 100 bp);
  - TP = truth events with at least one record in tp-base.vcf.gz; FN = truth events - TP;
  - FP = records of fp.vcf.gz, with mate pairs (MATEID/MATE_ID or reciprocal ALTs) counted once, and
    without the records whose mate is a TP call in tp-comp.vcf.gz (the second half of a correctly
    called event that truvari left unmatched).

Usage: truvari_events.py TRUTH.vcf TRUVARI_DIR      prints: events TP FP FN
"""
import re
import sys

import vcf


def read(path):
    return {r.ID: r for r in vcf.Reader(filename=path, compressed=path.endswith(".gz"))}


def mates_by_id(records):
    mates = {}
    for rid, r in records.items():
        for key in ("MATEID", "MATE_ID"):
            if key in r.INFO:
                m = r.INFO[key]
                m = m[0] if isinstance(m, (list, tuple)) else m
                mates[rid], mates[m] = m, rid
                break
    return mates


def alt_locus(r):
    for alt in r.ALT or []:
        chrom, pos = getattr(alt, "chr", None), getattr(alt, "pos", None)
        if chrom is not None and pos is not None:
            return str(chrom), int(pos)
        m = re.search(r"[\[\]]([^\[\]:]+):(\d+)[\[\]]", str(alt))
        if m:
            return m.group(1), int(m.group(2))
    return None


def mates_by_geometry(records, tol=100):
    """BND records whose ALTs name each other's locus, for callers that write no MATEID."""
    ends = [(rid, r.CHROM, int(r.POS)) + loc for rid, r in records.items() if (loc := alt_locus(r))]
    mates = {}
    for i, (a, ca, pa, ta, qa) in enumerate(ends):
        if a in mates:
            continue
        for b, cb, pb, tb, qb in ends[i + 1:]:
            if b in mates or ta != cb or tb != ca or abs(qa - pb) > tol or abs(qb - pa) > tol:
                continue
            mates[a], mates[b] = b, a
            break
    return mates


def mates(records):
    """MATEID pairs, completed by geometric pairs for the records that have none."""
    m = mates_by_id(records)
    for k, v in mates_by_geometry(records).items():
        m.setdefault(k, v)
    return m


def count_events(records, mate):
    counted, n = set(), 0
    for rid in records:
        if rid in counted:
            continue
        n += 1
        counted.add(rid)
        if mate.get(rid) in records:
            counted.add(mate[rid])
    return n


def main():
    truth_path, d = sys.argv[1:3]
    truth = read(truth_path)
    tp_base, fp, tp_comp = read(f"{d}/tp-base.vcf.gz"), read(f"{d}/fp.vcf.gz"), read(f"{d}/tp-comp.vcf.gz")

    truth_mates = mates(truth)
    events = count_events(truth, truth_mates)
    tp = count_events(tp_base, truth_mates)

    both = dict(fp)
    both.update(tp_comp)
    both_mates = mates(both)
    fp_kept = {rid: r for rid, r in fp.items() if both_mates.get(rid) not in tp_comp}
    n_fp = count_events(fp_kept, mates(fp))

    print(events, tp, n_fp, events - tp)


if __name__ == "__main__":
    main()
