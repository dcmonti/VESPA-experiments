#!/usr/bin/env python3
"""
Multi-caller benchmark PER EVENT against several GTs, with the event and match
definition of count_events.py (hg008/analysis_dcm/gates_2026-08-15/), which in turn
runs the one established in minda_ensemble/concordance_hifi/build_upset_hifi_eventlevel.py.

Why not bench_multicaller.sh: that goes through truvari + filter_mate_event.py,
which counts per set (tp-base, fn), collapsing the mates within each. A GT event
hit at one end only ends up ONCE AS TP AND ONCE AS FN. Here the event
is the breakpoint pair and a TP requires both ends within the window.

On top of the count_events definition, the two rules of the old bench are added,
which count_events does not have:

  1. length threshold MINLEN (default 50) on the EVENTS, identical for GT and calls:
     VESPA does not call below --min-sv-length, so without a threshold callers that
     emit 20-49 bp would be rewarded or penalized for the minimum cut. INS and
     intra-chromosomal pairs with known length < MINLEN are discarded; translocations
     and events of unknown length are kept (as truvari does with BNDs).
  2. two representation normalizations, both checked on all VCFs:
     - nanomonsv writes INS as <INS> with the length only in SVINSLEN: without
       SVLEN the length would be 0, the event would escape the threshold and the match
       would fall back to position only. SVLEN is derived from SVINSLEN.
     - Severus v1.0 writes INV with END=POS and the distance in SVLEN: without
       correction the two breakpoints of the event coincide. END becomes POS+|SVLEN|.
       The only caller with this convention; the others and the GT have END > POS.

The PASS filter (--caller-filter PASS) applies ONLY to the caller VCF, never to the
GT, which would otherwise be left with a different denominator.

USAGE:
  bench_events_multicaller.py --out BENCH.tsv [--minlen 50] \\
      --truth NAME=GT.vcf [--truth ...] --call NAME=VCF [--call ...]

Check: with --minlen 0 and caller-filter ALL it reproduces count_events.py.
"""
import argparse
import os
import re
import tempfile
from collections import defaultdict

A = os.environ["HG008"]
SRC = f"{A}/minda_ensemble/concordance_hifi/build_upset_hifi_eventlevel.py"
MARK = "# ---- costruisci eventi per tool ----"


def load_established():
    """parse_events/match from the established script, as count_events.py does: the
    file has no __main__ guard, so only the part before the marker is executed."""
    head = open(SRC).read().split(MARK)[0]
    ns = {"__name__": "_established"}
    exec(compile(head, SRC, "exec"), ns)
    return ns["parse_events"], ns["match"]


SVINSLEN = re.compile(r"(?:^|;)SVINSLEN=(\d+)")
HAS_SVLEN = re.compile(r"(?:^|;)SVLEN=")
SVLEN_RE = re.compile(r"(?:^|;)SVLEN=(-?\d+)")
END_RE = re.compile(r"(?:^|;)END=(\d+)")
SVTYPE_RE = re.compile(r"(?:^|;)SVTYPE=([A-Za-z:]+)")


def prepare(path, passonly):
    """Temporary copy of the VCF: PASS only if requested, plus the two
    representation normalizations. Returns
    (path, n_records, n_svlen_derived, n_end_fixed)."""
    fd, tmp = tempfile.mkstemp(suffix=".vcf")
    n = fixed = endfix = 0
    with os.fdopen(fd, "w") as out, open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                out.write(line)
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 8:
                continue
            if passonly and c[6] not in ("PASS", "."):
                continue
            if "SVTYPE=INS" in c[7] and not HAS_SVLEN.search(c[7]):
                m = SVINSLEN.search(c[7])
                if m:
                    c[7] = f"SVLEN={m.group(1)};{c[7]}"
                    fixed += 1
            # Severus v1.0 writes inversions with END=POS and the distance between the
            # breakpoints in SVLEN. parse_events uses END when present, so without
            # this correction ALL its INV become events with the two
            # breakpoints coinciding: they lose the true matches and, with a large SVLEN,
            # the window max(500, 0.3*span) produces false ones (a 135 Mb INV
            # ends up matching a 54 bp DEL 2 Mb away).
            # All the other callers and the GT write END > POS.
            # Only the types that really have a second breakpoint: for an INS
            # END=POS is the correct encoding (it consumes no reference) and
            # parse_events does not even look at it.
            msv = SVTYPE_RE.search(c[7])
            if msv and msv.group(1).upper() not in ("INS", "BND") \
                    and "[" not in c[4] and "]" not in c[4]:
                me, ml = END_RE.search(c[7]), SVLEN_RE.search(c[7])
                if me and ml and int(me.group(1)) == int(c[1]) and abs(int(ml.group(1))) > 0:
                    end = int(c[1]) + abs(int(ml.group(1)))
                    c[7] = c[7][:me.start(1)] + str(end) + c[7][me.end(1):]
                    endfix += 1
            n += 1
            out.write("\t".join(c) + "\n")
    return tmp, n, fixed, endfix


def size_ok(e, minlen):
    """Event within the threshold? Translocations and unknown lengths are kept."""
    if minlen <= 0 or not e["svlen"]:
        return True
    if e["kind"] == "PAIR" and e["bps"][0][0] != e["bps"][1][0]:
        return True
    return e["svlen"] >= minlen


def score(gt, cl, match):
    """Identical to the count of count_events.py."""
    def gkey(e):
        return ("INS", e["bps"][0][0]) if e["kind"] == "INS" else \
               ("PAIR", e["bps"][0][0], e["bps"][1][0])
    by = defaultdict(list)
    for i, e in enumerate(cl):
        by[gkey(e)].append(i)
    hit_gt = [False] * len(gt)
    hit_cl = [False] * len(cl)
    for gi, ge in enumerate(gt):
        for ci in by.get(gkey(ge), ()):
            if match(ge, cl[ci]):
                hit_gt[gi] = hit_cl[ci] = True
    tp = sum(hit_gt)
    return tp, len(cl) - sum(hit_cl), len(gt) - tp


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", action="append", required=True, help="NAME=GT.vcf")
    ap.add_argument("--call", action="append", required=True, help="NAME=VCF")
    ap.add_argument("--minlen", type=int, default=50)
    ap.add_argument("--caller-filter", choices=["ALL", "PASS", "BOTH"], default="BOTH")
    ap.add_argument("--regions-bed", default="",
                    help="BED of the GIAB benchmark regions: keeps only the events "
                         "with ALL anchors inside, both in the truth and "
                         "in the calls (equivalent of truvari --includebed). "
                         "Filters EVENTS, not records: filtering records would split "
                         "BND pairs.")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    parse_events, match = load_established()

    keep_region = None
    if a.regions_bed:
        from bisect import bisect_right
        reg = {}
        with open(a.regions_bed) as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track")):
                    continue
                c = line.split("\t")
                reg.setdefault(c[0], []).append((int(c[1]), int(c[2])))
        for k in reg:
            reg[k].sort()
        starts = {k: [x[0] for x in v] for k, v in reg.items()}

        def _inside(chrom, pos):
            v = reg.get(chrom)
            if not v:
                return False
            j = bisect_right(starts[chrom], pos - 1) - 1
            return j >= 0 and v[j][1] > pos - 1

        def keep_region(e):
            return all(_inside(c, p) for c, p in e["bps"])
    filters = ["ALL", "PASS"] if a.caller_filter == "BOTH" else [a.caller_filter]

    truths = {}
    for spec in a.truth:
        name, path = spec.split("=", 1)
        tmp, _, _, _ = prepare(path, passonly=False)
        truths[name] = [e for e in parse_events(tmp)
                        if size_ok(e, a.minlen)
                        and (keep_region is None or keep_region(e))]
        os.unlink(tmp)

    calls = {}
    for spec in a.call:
        name, path = spec.split("=", 1)
        for f in filters:
            tmp, n, fixed, endfix = prepare(path, passonly=(f == "PASS"))
            ev = [e for e in parse_events(tmp)
                  if size_ok(e, a.minlen)
                  and (keep_region is None or keep_region(e))]
            os.unlink(tmp)
            calls[(name, f)] = (ev, n, fixed, endfix)

    cols = ["truth", "truth_events", "caller_filter", "tool", "n_records",
            "n_svlen_from_svinslen", "n_end_fixed", "call_events", "TP", "FP",
            "FN", "precision", "recall", "f1"]
    with open(a.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for tname, gt in truths.items():
            for (cname, f), (ev, n, fixed, endfix) in calls.items():
                tp, fp, fn = score(gt, ev, match)
                p = tp / (tp + fp) if tp + fp else 0.0
                r = tp / (tp + fn) if tp + fn else 0.0
                f1 = 2 * p * r / (p + r) if p + r else 0.0
                out.write("\t".join(map(str, [
                    tname, len(gt), f, cname, n, fixed, endfix, len(ev), tp, fp, fn,
                    f"{p:.4f}", f"{r:.4f}", f"{f1:.4f}"])) + "\n")
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
