#!/usr/bin/env python3
"""HG008 evaluation rule: a call that matches a truth event outside the benchmark regions is neither a TP
nor an FP.

The calls of a Minda run are filtered with the region BED at the exact base, while Minda matches them
within its tolerance: a call a few bases outside an excluded gap is kept and, without this rule, counted
as an FP against a truth event the BED removed.

For a Minda run directory (callers.tsv, gt.vcf of the region, out/results/overall.tsv), Minda is run
again with the same calls against the whole truth set; the calls that match a truth event absent from
the region gt.vcf are removed from the FP count. TP and FN do not change.
Writes out/results/overall_rule.tsv: caller, TP, FN, FP, precision, recall, F1, IDs of the removed calls.

Usage: minda_excluded_truth.py RUN_DIR TRUTH_ALL.vcf MINDA_PYTHON MINDA_PY
"""
import csv
import os
import subprocess
import sys


def main():
    run, truth_all, python, minda = sys.argv[1:5]
    region_ids = {l.split("\t")[2] for l in open(f"{run}/gt.vcf") if not l.startswith("#")}
    od = f"{run}/out_alltruth"
    subprocess.run([python, minda, "truthset", "--base", truth_all, "--tsv", f"{run}/callers.tsv",
                    "--out_dir", od, "--tolerance", "500", "--min_size", "50", "--sample_name", "HG008"],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    rows = list(csv.DictReader(open(f"{run}/out/results/overall.tsv"), delimiter="\t"))
    callers = [r["Caller"] for r in rows if r["Caller"] != "base"]
    excluded = {c: [] for c in callers}
    for x in csv.DictReader(open(f"{od}/HG008_support.tsv"), delimiter="\t"):
        if x["ID_x"] in region_ids or x["ID_y"] in region_ids:
            continue
        for c in callers:
            if x.get(c) == "True":
                excluded[c].append(x["ID_x"])

    with open(f"{run}/out/results/overall_rule.tsv", "w") as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(["Caller", "True Positives", "False Negatives", "False Positives", "Precision", "Recall",
                    "F1 Score", "FP matched to excluded truth"])
        for r in rows:
            tp, fn = int(r["True Positives"]), int(r["False Negatives"])
            fp = int(r["False Positives"]) - len(excluded.get(r["Caller"], []))
            p = tp / (tp + fp) if tp + fp else 0
            rc = tp / (tp + fn) if tp + fn else 0
            f1 = 2 * p * rc / (p + rc) if p + rc else 0
            w.writerow([r["Caller"], tp, fn, fp, p, rc, f1, ",".join(excluded.get(r["Caller"], [])) or "."])


if __name__ == "__main__":
    main()
