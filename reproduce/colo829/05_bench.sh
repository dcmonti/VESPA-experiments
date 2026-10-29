#!/usr/bin/env bash
# COLO829: benchmark of the calls written by 04_call.sh, and of the published HiFi callsets of Qin & Li, against
# the COLO829 somatic SV truth set, with Minda and truvari.
#
# Truth: Valle-Inclan et al. 2022, lifted to GRCh38 (truthset_COLO829_hg38lifted.symbolic.vcf), genome-wide: there
# are no benchmark regions. Minimum SV size 50 bp.
# Records on contigs other than chr1-22, X, Y are dropped, in the truth and in the calls (bench/prep_vcf.py
# with a BED of the whole chromosomes). Calls: FILTER PASS only.
# Minda:   truthset, tolerance 500 bp, minimum size 50 bp.
# truvari: bench, refdist 500, pctsize 0.7, pctseq 0, any SV type (--typeignore), BND distance 100, minimum size
#          50 bp; the two records of a BND pair count as one event (bench/truvari_events.py).
# Published callsets: Qin & Li, Zenodo 10.5281/zenodo.17317401 (v0.2, COLO829_hifi1): the raw calls of Severus,
# SAVANA, nanomonsv and Sniffles2 (minimap2 -s50), and minisv (generation5_filterasm); SAVANA 1.3.6 filtered by
# minisv with the de novo assembly of the normal (l+s, threshold 3; Zenodo 10.5281/zenodo.21852833, v4).
# Keskus et al. 2025 (Severus paper), Zenodo 10.5281/zenodo.14541057 (v2.0): Severus 1.2 on COLO829 PacBio HiFi
# (original_vcfs/severus_pb.vcf.gz, no VAF filter). Both on the same Revio 2023Q2 reads.
#
# Output, in $OUT/bench/min<size>/:
#   summary.tsv   method, call set, TP, FP, FN, precision, recall, F1 (every run of every caller)
#   best.tsv      for each method and caller: the default run and the threshold with the highest F1
#
# Usage: bash 05_bench.sh
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"

OUT=${OUT:-$COLO829/reads_ec_align}
REF=$OUT/ref/GRCh38_GIABv3.chr1-22XY.fa
CALLS=$OUT/calls
BENCH=$OUT/bench
TRUTH=$COLO829/gt/truthset_COLO829_hg38lifted.symbolic.vcf
QINLI=$COLO829/external_qinli_2026-09-24/v02/COLO829_hifi1
QINLI4=$COLO829/external_qinli_v4_2026-10-08
KESKUS=$COLO829/external_severus_2026-10-08/extract/variant_calls_and_benchmarks/COLO829/original_vcfs

export PATH=$SV_BIN:$PATH                         # python (PyVCF), bcftools, bgzip, tabix, truvari 5.4.0
HELPERS=$(cd "$(dirname "$0")/../bench" && pwd)

THRESHOLDS="3 5 8 10 15 20 25"
MIN_SIZES="50"
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
  published QINLI-SEVERUS   $QINLI/severus/severus_somatic.vcf
  published QINLI-SAVANA    $QINLI/savana/grch38_T_tag.classified.somatic.vcf
  published QINLI-NANOMONSV $QINLI/nanomonsv/grch38_tnpair.vcf
  published QINLI-SNIFFLES  $QINLI/sniffles2/grch38_snf_somatic.vcf
  published QINLI-MINISV    $QINLI/minisv/minisv_generation5_filterasm.converted.vcf
  published QLMSV-SAVANA    $QINLI4/ensembl__COLO829_hifi1__savana_l+s_3_filtered.vcf
  published KESKUS-SEVERUS  $KESKUS/severus_pb.vcf.gz
}

# truvari bench of one call set; writes events, TP, FP, FN to <outdir>.counts.
# truvari merges the two headers, so the truth header first gets the INFO definitions of the calls.
truvari_one() { # min_size name
  local min=$1 name=$2 d=$BENCH/min$1
  local truth=$d/gt.vcf calls=$d/vcf/$name.vcf o=$d/truvari/$name
  rm -rf $o $o.tmp; mkdir -p $o.tmp
  { bcftools view -h $truth | grep -v '^#CHROM'; bcftools view -h $calls | grep '^##INFO'; bcftools view -h $truth | grep '^#CHROM'; } > $o.tmp/header
  bcftools view $truth | bgzip -c > $o.tmp/truth0.vcf.gz
  bcftools reheader -h $o.tmp/header -o $o.tmp/truth.vcf.gz $o.tmp/truth0.vcf.gz
  bgzip -c $calls > $o.tmp/calls.vcf.gz
  tabix -p vcf $o.tmp/truth.vcf.gz; tabix -p vcf $o.tmp/calls.vcf.gz
  truvari bench -b $o.tmp/truth.vcf.gz -c $o.tmp/calls.vcf.gz -f $REF -o $o \
    --sizemin $min --sizefilt $min --sizemax 100000000 --pctseq 0 --typeignore > $o.log 2>&1
  rm -rf $o.tmp
  read -r _ tp fp fn < <(python $HELPERS/truvari_events.py $truth $o)
  echo -e "$name\t$tp\t$fp\t$fn" > $o.counts
}
export -f truvari_one
export BENCH REF HELPERS

mkdir -p $BENCH
awk -v OFS='\t' '{print $1, 0, $2}' $REF.fai > $BENCH/chromosomes.bed

for min in $MIN_SIZES; do
  d=$BENCH/min$min
  rm -rf $d; mkdir -p $d/vcf $d/truvari
  python $HELPERS/prep_vcf.py $TRUTH $d/gt.vcf --bed $BENCH/chromosomes.bed --min-len $min

  # calls: PASS only, on chr1-22, X, Y
  : > $d/callers.tsv
  while read -r name vcf; do
    python $HELPERS/prep_vcf.py <(zcat -f $vcf) $d/vcf/$name.vcf --bed $BENCH/chromosomes.bed --pass-only
    echo -e "$d/vcf/$name.vcf\t$name" >> $d/callers.tsv
  done < <(callsets)

  # Minda
  $MINDA_PYTHON $MINDA truthset --base $d/gt.vcf --tsv $d/callers.tsv --out_dir $d/out \
    --tolerance 500 --min_size $min --sample_name COLO829 > $d/minda.log 2>&1

  # truvari, PARALLEL call sets at a time
  cut -f2 $d/callers.tsv | xargs -P $PARALLEL -I{} bash -c "set -euo pipefail; truvari_one $min {}"

  # summary
  {
    echo -e "method\tcallset\tTP\tFP\tFN\tprecision\trecall\tF1"
    awk -F'\t' 'NR > 1 {print "minda\t" $1 "\t" $2 "\t" $4 "\t" $3}' $d/out/results/overall.tsv
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
