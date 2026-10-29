#!/usr/bin/env bash
# Germline-removal cascade (Results, block A): Minda on the call sets of each row, both region sets,
# with the excluded-truth rule, exactly as reproduce/hg008/05_bench.sh (tolerance 500, min size 50).
# Rows are given as NAME VCF PASSONLY(0|1) in rows.tsv next to this script; a row whose VCF does not
# exist yet is skipped.
# Output: bench/<regions>/out/results/overall_rule.tsv and bench/summary.tsv (row, regions, TP, FP, FN, P, R, F1).
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"
D=$HG008/reads_ec_align/germline_cascade
# one Minda run at a time: 03_call_tables.sh and 06_v63.sh may both call this script
exec 9> $D/bench/.lock; flock 9
R=$HG008/reads_ec_align
V05=$HG008/ground_truths/v0.5
export PATH=$SV_BIN:$PATH
HELPERS=$REPRO/bench
for regions in clonal clonal_and_subclonal; do
  d=$D/bench/$regions; rm -rf $d; mkdir -p $d/vcf
  cp $R/bench/$regions/gt.vcf $d/gt.vcf           # same truth as the caller comparison
  : > $d/callers.tsv
  while read -r name vcf passonly; do
    [[ -z $name || $name == \#* ]] && continue
    [[ -s $vcf ]] || { echo "skip $name: $vcf missing" >&2; continue; }
    flag=""; [[ $passonly == 1 ]] && flag="--pass-only"
    python $HELPERS/prep_vcf.py $vcf $d/vcf/$name.vcf --bed $V05/GRCh38_HG008-T-V0.5_somatic-stvar-$regions.draftbenchmark.bed $flag
    echo -e "$d/vcf/$name.vcf\t$name" >> $d/callers.tsv
  done < $D/rows.tsv
  $MINDA_PYTHON $MINDA truthset --base $d/gt.vcf --tsv $d/callers.tsv --out_dir $d/out \
    --tolerance 500 --min_size 50 --sample_name HG008 > $d/minda.log 2>&1
  python $HELPERS/minda_excluded_truth.py $d $R/bench/truth_all.vcf $MINDA_PYTHON $MINDA
done
{ echo -e "row\tregions\tTP\tFP\tFN\tprecision\trecall\tF1"
  for regions in clonal clonal_and_subclonal; do
    awk -F'\t' -v r=$regions 'NR>1 && $1!="base" {printf "%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\n", $1, r, $2, $4, $3, $5, $6, $7}' $D/bench/$regions/out/results/overall_rule.tsv
  done; } > $D/bench/summary.tsv
cat $D/bench/summary.tsv
