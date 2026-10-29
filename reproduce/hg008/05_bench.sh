#!/usr/bin/env bash
# HG008: benchmark of the calls written by 04_call.sh, and of the published HiFi callsets of the Minda ensemble,
# against the GIAB HG008-T V0.5 draft somatic SV benchmark, with Minda and truvari.
#
# Truth: GRCh38_HG008-T-V0.5_somatic-stvar_ALL.draftbenchmark.vcf, SVs of at least 50 bp, in two region sets
# (clonal, clonal_and_subclonal). Calls: FILTER PASS only. In both, a record is kept only if all its anchors
# fall in the regions (bench/prep_vcf.py).
# Minda:   truthset, tolerance 500 bp, minimum size 50 bp.
# truvari: bench, refdist 500, pctsize 0.7, pctseq 0, any SV type (--typeignore), BND distance 100, SVs of at
#          least 50 bp; the two records of a BND pair count as one event (bench/truvari_events.py).
# Rule for both: a call that matches a truth event outside the regions is neither a TP nor an FP
# (Minda: bench/minda_excluded_truth.py; truvari: FP from a second run against the whole truth set).
#
# Output, in $OUT/bench/<regions>/:
#   summary.tsv   method, call set, TP, FP, FN, precision, recall, F1 (every run of every caller)
#   best.tsv      for each method and caller: the default run and the threshold with the highest F1
#
# Usage: bash 05_bench.sh
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"

OUT=${OUT:-$HG008/reads_ec_align}
REF=$OUT/ref/GRCh38_GIABv3.chr1-22X.fa
CALLS=$OUT/calls
BENCH=$OUT/bench
V05=$HG008/ground_truths/v0.5
TRUTH=$V05/GRCh38_HG008-T-V0.5_somatic-stvar_ALL.draftbenchmark.vcf
# Published HiFi callsets of the Minda ensemble (NIH-NCI_minda-ensemble_20240710/vcf)
PUBLISHED=$HG008/minda_ensemble/vcf

export PATH=$SV_BIN:$PATH                         # python (PyVCF), bcftools, bgzip, tabix, truvari 5.4.0
HELPERS=$(cd "$(dirname "$0")/../bench" && pwd)

THRESHOLDS="3 5 8 10 15 20 25"
PARALLEL=8

# Call sets: name and VCF. Minda attributes matches by name prefix, so no name is a prefix of another.
# published NAME VCF: a published call set, listed only if it has been downloaded (README).
published() { [[ -s $2 ]] && echo "$1 $2" || echo "skipped $1: $2 not found" >&2; }

callsets() {
  for run in default $(for t in $THRESHOLDS; do printf "t%02d " $t; done); do
    echo "SEVERUS_$run   $CALLS/severus/$run/somatic_SVs/severus_somatic.vcf"
    echo "SAVANA_$run    $CALLS/savana/$run/tumor.raw.hp.classified.somatic.vcf"
    echo "NANOMONSV_$run $CALLS/nanomonsv/$run/tumor.nanomonsv.result.vcf"
    echo "SNIFFLES_$run  $CALLS/sniffles2/$run/merge_somatic_only.vcf"
    [[ $run != default ]] && echo "VESPA_$run     $CALLS/vespa/$run/calls.vcf"
  done
  published PUBLISHED-SEVERUS   $PUBLISHED/severus_hifi_vaf.vcf
  published PUBLISHED-SAVANA    $PUBLISHED/savana_hifi_vaf.vcf
  published PUBLISHED-NANOMONSV $PUBLISHED/nanomonsv_hifi_vaf.vcf
  published PUBLISHED-SNIFFLES  $PUBLISHED/sniffles_hifi_vaf.vcf
}

# truvari bench of one call set against one truth VCF; prints events TP FP FN.
# truvari merges the two headers, so the truth header first gets the INFO definitions of the calls.
truvari_events() { # truth calls outdir
  local truth=$1 calls=$2 o=$3
  rm -rf $o $o.tmp; mkdir -p $o.tmp
  { bcftools view -h $truth | grep -v '^#CHROM'; bcftools view -h $calls | grep '^##INFO'; bcftools view -h $truth | grep '^#CHROM'; } > $o.tmp/header
  bcftools view $truth | bgzip -c > $o.tmp/truth0.vcf.gz
  bcftools reheader -h $o.tmp/header -o $o.tmp/truth.vcf.gz $o.tmp/truth0.vcf.gz
  bgzip -c $calls > $o.tmp/calls.vcf.gz
  tabix -p vcf $o.tmp/truth.vcf.gz; tabix -p vcf $o.tmp/calls.vcf.gz
  truvari bench -b $o.tmp/truth.vcf.gz -c $o.tmp/calls.vcf.gz -f $REF -o $o \
    --sizemin 50 --sizefilt 50 --sizemax 100000000 --pctseq 0 --typeignore > $o.log 2>&1
  rm -rf $o.tmp
  python $HELPERS/truvari_events.py $truth $o
}

# One call set in one region set: truvari against the region truth (TP, FN) and against the whole truth (FP).
truvari_one() { # regions name
  local d=$BENCH/$1 name=$2
  read -r _ tp _ fn < <(truvari_events $d/gt.vcf $d/vcf/$name.vcf $d/truvari/$name)
  read -r _ _ fp _  < <(truvari_events $BENCH/truth_all.vcf $d/vcf/$name.vcf $d/truvari_alltruth/$name)
  echo -e "$name\t$tp\t$fp\t$fn" > $d/truvari/$name.counts
}
export -f truvari_events truvari_one
export BENCH REF HELPERS

mkdir -p $BENCH
python $HELPERS/prep_vcf.py $TRUTH $BENCH/truth_all.vcf --min-len 50

for regions in clonal clonal_and_subclonal; do
  d=$BENCH/$regions
  rm -rf $d; mkdir -p $d/vcf $d/truvari $d/truvari_alltruth
  python $HELPERS/prep_vcf.py $TRUTH $d/gt.vcf --bed $V05/GRCh38_HG008-T-V0.5_somatic-stvar-$regions.draftbenchmark.bed --min-len 50

  # calls: PASS only, in the regions
  : > $d/callers.tsv
  while read -r name vcf; do
    python $HELPERS/prep_vcf.py $vcf $d/vcf/$name.vcf --bed $V05/GRCh38_HG008-T-V0.5_somatic-stvar-$regions.draftbenchmark.bed --pass-only
    echo -e "$d/vcf/$name.vcf\t$name" >> $d/callers.tsv
  done < <(callsets)

  # Minda, then the excluded-truth rule
  $MINDA_PYTHON $MINDA truthset --base $d/gt.vcf --tsv $d/callers.tsv --out_dir $d/out \
    --tolerance 500 --min_size 50 --sample_name HG008 > $d/minda.log 2>&1
  python $HELPERS/minda_excluded_truth.py $d $BENCH/truth_all.vcf $MINDA_PYTHON $MINDA

  # truvari, PARALLEL call sets at a time
  cut -f2 $d/callers.tsv | xargs -P $PARALLEL -I{} bash -c "set -euo pipefail; truvari_one $regions {}"

  # summary
  {
    echo -e "method\tcallset\tTP\tFP\tFN\tprecision\trecall\tF1"
    awk -F'\t' 'NR > 1 {print "minda\t" $1 "\t" $2 "\t" $4 "\t" $3}' $d/out/results/overall_rule.tsv
    cat $d/truvari/*.counts | awk -F'\t' '{print "truvari\t" $1 "\t" $2 "\t" $3 "\t" $4}'
  } | awk -F'\t' -v OFS='\t' 'NR == 1 {print; next}
        { p = $3 + $4 ? $3 / ($3 + $4) : 0; r = $3 + $5 ? $3 / ($3 + $5) : 0
          printf "%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\n", $1, $2, $3, $4, $5, p, r, p + r ? 2 * p * r / (p + r) : 0 }' \
    > $d/summary.tsv

  # best threshold per method and caller (highest F1 over t03..t25), and the default run
  { head -1 $d/summary.tsv
    awk -F'\t' 'NR > 1 { split($2, a, "_"); key = $1 "\t" a[1] }
          NR > 1 && (a[2] == "default" || $2 !~ /_/) {print; next}
          NR > 1 && (!(key in f) || $8 > f[key]) {f[key] = $8; best[key] = $0}
          END {for (k in best) print best[k]}' $d/summary.tsv | sort -k1,1 -k2,2
  } > $d/best.tsv
done
