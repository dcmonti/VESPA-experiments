#!/usr/bin/env bash
# Effect of each germline filter, HG008 (w = 5) and COLO829 -s50 (w = 8): the runs of the main Results.
# Rows, cumulative in the order the filters are applied (Methods, Assessing the somatic origin of candidates):
#   F0_none         VESPA call without the normal: no germline filter
#   F1_sameallele   the call with the normal, every record: only the same-allele test, which removes
#   F2_germline     ... without the records flagged GERMLINE
#   F3_polynormal   ... and POLYNORMAL
#   F4_repeatdepth  ... and REPEATDEPTH
#   F5_lowpower     ... and LOWPOWER: the PASS records, = calls/vespa/tNN of the caller comparison
# Each record carries one flag, the first filter it fails, so a row removes what that filter flags after the
# ones before it.
# Calls: HG008 reuses reads_ec_align/bp_coverage/{nogl_w5,with_normal_w5} (flags of 04_call.sh, w = 5; with the
# normal = calls/vespa/t05). COLO829 runs one VESPA call without the normal, flags of 04_call.sh, w = 8, on the
# -s50 alignments; with the normal = reads_ec_align_s50/calls/vespa/t08.
# Bench: Minda as 05_bench.sh of each dataset (tolerance 500, min size 50); HG008 within the clonal and the
# clonal and subclonal regions of V0.5 with the excluded-truth rule, COLO829 genome-wide on chr1-22, X, Y.
# Output: <dataset run>/germline_filters/{rows,bench,summary.tsv}.  Usage: bash run.sh [hg008] [colo829]
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"
export PATH=$SV_BIN:$PATH
HERE=$(cd "$(dirname "$0")" && pwd)
HELPERS=$(cd $HERE/../bench && pwd)
EXCLUDE=$REFS/grch38_centromers.bed
FLAGS="GERMLINE POLYNORMAL REPEATDEPTH LOWPOWER"
NAMES="F2_germline F3_polynormal F4_repeatdepth F5_lowpower"

# rows NOFILTERS_VCF WITHNORMAL_VCF OUTDIR: one VCF per row in OUTDIR/rows, list in OUTDIR/rows.tsv
rows() {
  local none=$1 full=$2 o=$3 drop="" i=0
  mkdir -p $o/rows
  cp $none $o/rows/F0_none.vcf; cp $full $o/rows/F1_sameallele.vcf
  set -- $NAMES
  for f in $FLAGS; do
    drop="$drop $f"
    awk -F'\t' -v drop="$drop" 'BEGIN {n = split(drop, d, " "); for (i = 1; i <= n; i++) x[d[i]] = 1}
      /^#/ || !($7 in x)' $full > $o/rows/$1.vcf
    shift
  done
  for n in F0_none F1_sameallele $NAMES; do echo -e "$n\t$o/rows/$n.vcf"; done > $o/rows.tsv
}

hg008() {
  local R=$HG008/reads_ec_align
  local V05=$HG008/ground_truths/v0.5
  local o=$R/germline_filters
  rows $R/bp_coverage/nogl_w5/calls.vcf $R/bp_coverage/with_normal_w5/calls.vcf $o
  for regions in clonal clonal_and_subclonal; do
    local d=$o/bench/$regions; rm -rf $d; mkdir -p $d/vcf
    cp $R/bench/$regions/gt.vcf $d/gt.vcf; : > $d/callers.tsv
    while read -r name vcf; do
      python $HELPERS/prep_vcf.py $vcf $d/vcf/$name.vcf --bed $V05/GRCh38_HG008-T-V0.5_somatic-stvar-$regions.draftbenchmark.bed
      echo -e "$d/vcf/$name.vcf\t$name" >> $d/callers.tsv
    done < $o/rows.tsv
    $MINDA_PYTHON $MINDA truthset --base $d/gt.vcf --tsv $d/callers.tsv --out_dir $d/out --tolerance 500 --min_size 50 \
      --sample_name HG008 > $d/minda.log 2>&1
    python $HELPERS/minda_excluded_truth.py $d $R/bench/truth_all.vcf $MINDA_PYTHON $MINDA
  done
  { echo -e "row\tregions\tTP\tFP\tFN\tprecision\trecall\tF1"
    for regions in clonal clonal_and_subclonal; do
      awk -F'\t' -v r=$regions 'NR > 1 && $1 != "base" {printf "%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\n", $1, r, $2, $4, $3, $5, $6, $7}' \
        $o/bench/$regions/out/results/overall_rule.tsv
    done; } > $o/summary.tsv
  cat $o/summary.tsv
}

colo829() {
  local R=$COLO829/reads_ec_align_s50
  local o=$R/germline_filters c=$R/germline_filters/calls/nofilters_t08
  local REF=$R/ref/GRCh38_GIABv3.chr1-22XY.fa
  local TRUTH=$COLO829/gt/truthset_COLO829_hg38lifted.symbolic.vcf
  mkdir -p $c
  if [[ ! -s $c/calls.vcf ]]; then
    $VESPA call --reference $REF --bam $R/align/tumor.ec.hp.bam --sfs $R/search/tumor.dedup.cap500.sfs --bed-exclusion $EXCLUDE \
      --min-cluster-weight 8 --min-sv-length 50 --min-mapq 60 --clipped --no-clipped-fallback --threads 32 \
      > $c/calls.unsorted.vcf 2> $c/call.log
    bcftools sort $c/calls.unsorted.vcf -o $c/calls.vcf && rm $c/calls.unsorted.vcf
  fi
  rows $c/calls.vcf $R/calls/vespa/t08/calls.vcf $o
  local d=$o/bench; rm -rf $d; mkdir -p $d/vcf
  python $HELPERS/prep_vcf.py $TRUTH $d/gt.vcf --bed $R/bench/chromosomes.bed --min-len 50
  : > $d/callers.tsv
  while read -r name vcf; do
    python $HELPERS/prep_vcf.py $vcf $d/vcf/$name.vcf --bed $R/bench/chromosomes.bed
    echo -e "$d/vcf/$name.vcf\t$name" >> $d/callers.tsv
  done < $o/rows.tsv
  $MINDA_PYTHON $MINDA truthset --base $d/gt.vcf --tsv $d/callers.tsv --out_dir $d/out --tolerance 500 --min_size 50 \
    --sample_name COLO829 > $d/minda.log 2>&1
  { echo -e "row\tTP\tFP\tFN\tprecision\trecall\tF1"
    awk -F'\t' 'NR > 1 {tp = $2; fn = $3; fp = $4; p = tp + fp ? tp / (tp + fp) : 0; r = tp + fn ? tp / (tp + fn) : 0
      printf "%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\n", $1, tp, fp, fn, p, r, p + r ? 2 * p * r / (p + r) : 0}' $d/out/results/overall.tsv
  } > $o/summary.tsv
  cat $o/summary.tsv
}

for ds in ${@:-hg008 colo829}; do $ds; done
