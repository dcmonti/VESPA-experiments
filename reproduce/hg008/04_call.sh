#!/usr/bin/env bash
# HG008: somatic SV calling with VESPA, Severus, SAVANA, nanomonsv and Sniffles2.
#
# Severus, SAVANA, nanomonsv, Sniffles2: the commands of the GIAB HG008 Minda ensemble
# (NIH-NCI_minda-ensemble_20240710, individual-variant-caller-runs/scripts/*_HG008_HiFi*.sh), run with the
# current release of each tool:
#   default   the Minda command
#   tNN       the Minda command with the tool's minimum number of tumor supporting reads set to NN
# VESPA: the same flags for every run, with the minimum cluster weight set to NN (tNN).
# The four other callers read the raw-read alignments, VESPA the error-corrected alignments and the SFS.
#
# Every command is timed with GNU time (/usr/bin/time -v, reports in $OUT/time/). The default runs, and VESPA
# at t05 (the minimum cluster weight of the pipeline configuration), run one at a time, so that their times
# are not affected by other jobs, with 32 threads each; then the grid runs, the five callers in parallel with
# 16 threads each.
#
# Output, in $OUT/calls/<caller>/<run>/:
#   severus     somatic_SVs/severus_somatic.vcf
#   savana      tumor.raw.hp.classified.somatic.vcf
#   nanomonsv   tumor.nanomonsv.result.vcf
#   sniffles2   merge_somatic_only.vcf
#   vespa       calls.vcf
#
# Usage: bash 04_call.sh
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"

OUT=${OUT:-$HG008/reads_ec_align}
REF=$OUT/ref/GRCh38_GIABv3.chr1-22X.fa
T_RAW=$OUT/align/tumor.raw.hp.bam
N_RAW=$OUT/align/normal.raw.hp.bam
T_EC=$OUT/align/tumor.ec.hp.bam
N_EC=$OUT/align/normal.ec.hp.bam
PHASED=$OUT/haplotag/normal.phased.vcf.gz
SFS=$OUT/search/tumor.dedup.cap500.sfs
# Tandem repeats: vntrs/human_GRCh38_no_alt_analysis_set.trf.bed of the Severus repository, the file used by Minda.
VNTR=$REFS/human_GRCh38_no_alt_analysis_set.trf.bed
# Centromeres and telomeres of GRCh38, excluded by VESPA.
EXCLUDE=$REFS/grch38_centromers.bed

THRESHOLDS="3 5 8 10 15 20 25"
THREADS_DEFAULT=32   # every caller, default runs (one at a time)
THREADS_GRID=16      # every caller, grid runs (five callers at a time)

mkdir -p $OUT/calls $OUT/logs $OUT/time

# timed NAME COMMAND...: run COMMAND under GNU time, report in $OUT/time/NAME.txt.
timed() { local name=$1; shift; /usr/bin/time -v -o $OUT/time/$name.txt "$@"; }
grid() { for t in $THRESHOLDS; do printf "t%02d\n" $t; done; }
threshold() { echo $((10#${1#t})); }

# Each caller function runs the runs given as arguments ("default" or tNN), one after the other.

severus() {
  for run in "$@"; do
    o=$OUT/calls/severus/$run
    [[ -s $o/somatic_SVs/severus_somatic.vcf ]] && continue
    extra=""; [[ $run != default ]] && extra="--min-support $(threshold $run)"
    rm -rf $o; mkdir -p $o
    timed severus.$run $SEVERUS --target-bam $T_RAW --control-bam $N_RAW --vntr-bed $VNTR --out-dir $o -t $THREADS \
      --phasing-vcf $PHASED --max-genomic-len 2000000 --TIN-ratio 0 $extra > $o/severus.log 2>&1
  done
}

savana() {
  for run in "$@"; do
    o=$OUT/calls/savana/$run
    [[ -s $o/tumor.raw.hp.classified.somatic.vcf ]] && continue
    extra=""; [[ $run != default ]] && extra="--min_support $(threshold $run)"
    rm -rf $o; mkdir -p $o
    timed savana.$run $SAVANA --tumour $T_RAW --normal $N_RAW --ref $REF --length 50 --threads $THREADS --outdir $o $extra > $o/savana.log 2>&1
  done
}

nanomonsv() {
  p=$OUT/calls/nanomonsv/parse
  mkdir -p $p
  [[ -s $p/tumor.bp_info.sorted.bed.gz ]]  || timed nanomonsv_parse.tumor  $NANOMONSV parse $T_RAW $p/tumor  > $p/parse_tumor.log 2>&1
  [[ -s $p/normal.bp_info.sorted.bed.gz ]] || timed nanomonsv_parse.normal $NANOMONSV parse $N_RAW $p/normal > $p/parse_normal.log 2>&1
  for run in "$@"; do
    o=$OUT/calls/nanomonsv/$run
    [[ -s $o/tumor.nanomonsv.result.vcf ]] && continue
    extra=""; [[ $run != default ]] && extra="--min_tumor_variant_read_num $(threshold $run)"
    rm -rf $o; mkdir -p $o
    # get writes its results next to the tumor prefix: each run gets its own links to the parse output
    ln -s $p/tumor.* $p/normal.* $o/
    timed nanomonsv_get.$run $NANOMONSV get $o/tumor $T_RAW $REF --control_prefix $o/normal --control_bam $N_RAW \
      --single_bnd --min_indel_size 50 --processes $THREADS $extra > $o/get.log 2>&1
    find $o -maxdepth 1 -type l -delete
  done
}

sniffles2() {
  s=$OUT/calls/sniffles2
  mkdir -p $s
  [[ -s $s/normal.snf ]] || timed sniffles.normal $SNIFFLES --input $N_RAW --snf $s/normal.snf --tandem-repeats $VNTR --reference $REF \
    --threads $THREADS --minsvlen 50 > $s/normal.log 2>&1
  for run in "$@"; do
    o=$s/$run
    [[ -s $o/merge_somatic_only.vcf ]] && continue
    extra=""; [[ $run != default ]] && extra="--minsupport $(threshold $run)"
    rm -rf $o; mkdir -p $o
    timed sniffles.tumor.$run $SNIFFLES --input $T_RAW --snf $o/tumor.snf --tandem-repeats $VNTR --reference $REF \
      --threads $THREADS --minsvlen 50 $extra > $o/tumor.log 2>&1
    # normal first, tumor second: SUPP_VEC '01' = called in the tumor, not in the normal
    timed sniffles_combine.$run $SNIFFLES --input $s/normal.snf $o/tumor.snf --vcf $o/merge_normal_tumor.vcf --threads $THREADS > $o/combine.log 2>&1
    timed sniffles_somatic.$run $BCFTOOLS view -i "SUPP_VEC = '01'" $o/merge_normal_tumor.vcf > $o/merge_somatic_only.vcf
  done
}

vespa() {
  for run in "$@"; do
    o=$OUT/calls/vespa/$run
    [[ -s $o/calls.vcf ]] && continue
    rm -rf $o; mkdir -p $o
    timed vespa.$run $VESPA call --reference $REF --bam $T_EC --sfs $SFS --normal-contigs-bam $N_EC --bed-exclusion $EXCLUDE \
      --min-cluster-weight $(threshold $run) --min-sv-length 50 --min-mapq 60 --clipped --no-clipped-fallback --threads $THREADS \
      > $o/calls.unsorted.vcf 2> $o/call.log
    timed vespa_sort.$run $BCFTOOLS sort $o/calls.unsorted.vcf -o $o/calls.vcf
    rm $o/calls.unsorted.vcf
  done
}

# 1. Default runs, one at a time.
THREADS=$THREADS_DEFAULT
severus default   > $OUT/logs/calls.severus.default.log   2>&1
savana default    > $OUT/logs/calls.savana.default.log    2>&1
nanomonsv default > $OUT/logs/calls.nanomonsv.default.log 2>&1
sniffles2 default > $OUT/logs/calls.sniffles2.default.log 2>&1
vespa t05         > $OUT/logs/calls.vespa.t05.log         2>&1

# 2. Grid: the five callers in parallel, the runs of each caller one after the other.
THREADS=$THREADS_GRID
for caller in severus savana nanomonsv sniffles2 vespa; do
  $caller $(grid) > $OUT/logs/calls.$caller.grid.log 2>&1 &
done
failed=0
for pid in $(jobs -p); do wait $pid || failed=1; done
[[ $failed == 0 ]] || { echo "a caller failed, see $OUT/logs/calls.*.log and the run logs" >&2; exit 1; }
