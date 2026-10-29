#!/usr/bin/env python3
"""Per-event VAF from our svviz2 runs (2026-10-27), joined with the Minda support matrix of the bench.

  HG008    svviz2_vaf/reports (V0.5 events of the clonal and subclonal benchmark, tumor raw reads)
           bench: bench/clonal_and_subclonal/out/HG008_support.tsv
  COLO829  colo829/reads_ec_align_s50/svviz2_vaf/reports (truth set, tumor raw reads aligned with -s50)
           bench: colo829/reads_ec_align_s50/bench/min50/out/COLO829_support.tsv

A report is named after the svviz2 variant: the record ID, or for a BND pair the ID of one mate (HG008) or the EVENT
(COLO829). ALT and REF are the svviz2 'count' rows. Events without a report get empty alt/ref/vaf.
walt, wref, wvaf: the svviz2 'weighted_count' rows truncated to integers, and wvaf = walt / (walt + wref). This is the
GIAB computation (SVVIZBYDATASET of V0.5 = ALT/total from floor(weighted_count)): on HG008 it reproduces the HiFi
totals within 2 reads for every event except two insertions (SV_219, SV_242), while 'count' overstates REF of every
tandem duplication by ~1.45x (checked 2026-10-27; whole-genome realignment gives the same counts as local).
Output: <dataset svviz2_vaf>/vaf_events.tsv with event, svtype, alt, ref, vaf, walt, wref, wvaf, giab_alt, giab_total
(HG008 only), subclonal, and the support columns of the bench (True/False per call set).
Usage: python svviz2_vaf_table.py
"""
import csv, glob, os, re
H, C = os.environ["HG008"], os.environ["COLO829"]
SETS = {"HG008": (f"{H}/reads_ec_align/svviz2_vaf",
                  f"{H}/reads_ec_align/bench/clonal_and_subclonal/out/HG008_support.tsv"),
        "COLO829": (f"{C}/reads_ec_align_s50/svviz2_vaf",
                    f"{C}/reads_ec_align_s50/bench/min50/out/COLO829_support.tsv")}


def counts(path):
    c = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["key"] in ("count", "weighted_count") and r["allele"] in ("alt", "ref"):
            c[r["key"], r["allele"]] = int(float(r["value"]))
    return tuple(c.get((k, a)) for k in ("count", "weighted_count") for a in ("alt", "ref"))


for name, (vdir, support) in SETS.items():
    rep = {os.path.basename(p).split(".")[0]: counts(p) for p in glob.glob(f"{vdir}/reports/*.report.tsv")}
    rows = [r for r in csv.DictReader(open(support), delimiter="\t") if r["base"] == "True"]
    runs = [k for k in rows[0] if k.split("_")[0] in ("VESPA", "SEVERUS", "SAVANA", "NANOMONSV", "SNIFFLES") or k.startswith(("PUBLISHED", "QINLI", "QLMSV", "KESKUS"))]
    out = []
    for r in rows:
        keys = [r["ID_x"], r["ID_y"]]
        for k in ("MATEID", "EVENT"):
            m = re.search(rf"(?:^|;){k}=([^;]*)", r["INFO_x"])
            if m: keys.append(m.group(1))
        alt, ref, walt, wref = next((rep[k] for k in keys if k in rep), (None,) * 4)
        vaf = alt / (alt + ref) if alt is not None and alt + ref > 0 else None
        wvaf = walt / (walt + wref) if walt is not None and walt + wref > 0 else None
        g = re.search(r"HG8T-PB-Hifi-116x:(\d+)/(\d+)", r["INFO_x"])
        out.append({"event": r["ID_x"], "svtype": r["SVTYPE"], "alt": "" if alt is None else alt,
                    "ref": "" if ref is None else ref, "vaf": "" if vaf is None else f"{vaf:.4f}",
                    "walt": "" if walt is None else walt, "wref": "" if wref is None else wref,
                    "wvaf": "" if wvaf is None else f"{wvaf:.4f}",
                    "giab_alt": g.group(1) if g else "", "giab_total": g.group(2) if g else "",
                    "subclonal": (re.search(r"SUBCLONAL=([^;]*)", r["INFO_x"]) or [None, ""])[1],
                    **{k: r[k] for k in runs}})
    with open(f"{vdir}/vaf_events.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t"); w.writeheader(); w.writerows(out)
    n = sum(o["vaf"] != "" for o in out)
    print(f"{name}: {len(out)} events, {n} with svviz2 VAF, {len(rep)} reports -> {vdir}/vaf_events.tsv")
