#!/usr/bin/env bash
# HG008, fourth arm of the 2x2 (homopolymer compression x SFS length cap): no compression AND no length cap.
# Uses the uncapped tumor SFS of run.sh (search/tumor.sfs, kept as search/tumor.sfs.keep, a hard link made while
# run.sh was running, because run.sh deletes tumor.sfs at the end), the index of run.sh, and the same call and bench.
#   1. dedup_sfs.py without a cap (search/tumor.dedup.nocap.sfs)
#   2. VESPA call at s = 5 (the default run of 04_call.sh), 32 threads -> calls_nocap/t05
#   3. Minda bench as 05_bench.sh (PASS, regions, excluded-truth rule) -> bench_nocap/<regions>
# Usage: bash run_nocap.sh
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"
R=$HG008/reads_ec_align
O=$R/nohpc
HERE=$(cd "$(dirname "$0")" && pwd)
HELPERS=$(cd $HERE/../bench && pwd)
V05=$HG008/ground_truths/v0.5
REF=$R/ref/GRCh38_GIABv3.chr1-22X.fa
EXCLUDE=$REFS/grch38_centromers.bed
timed() { local name=$1; shift; /usr/bin/time -v -o $O/time/$name.txt "$@"; }
export PATH=$SV_BIN:$PATH
[[ -s $O/search/tumor.sfs.keep ]] || { echo "missing $O/search/tumor.sfs.keep" >&2; exit 1; }
[[ -s $O/search/tumor.dedup.nocap.sfs ]] || timed dedup_nocap $PYTHON $LIB/dedup_sfs.py $O/search/tumor.sfs.keep \
  $O/search/tumor.dedup.nocap.sfs --max-length 0 2> $O/logs/dedup_nocap.log
o=$O/calls_nocap/t05
if [[ ! -s $o/calls.vcf ]]; then
  rm -rf $o; mkdir -p $o
  timed vespa_nocap.t05 $VESPA call --reference $REF --bam $R/align/tumor.ec.hp.bam --sfs $O/search/tumor.dedup.nocap.sfs \
    --normal-contigs-bam $R/align/normal.ec.hp.bam --bed-exclusion $EXCLUDE --min-cluster-weight 5 --min-sv-length 50 \
    --min-mapq 60 --clipped --no-clipped-fallback --threads 32 > $o/calls.unsorted.vcf 2> $o/call.log
  $BCFTOOLS sort $o/calls.unsorted.vcf -o $o/calls.vcf && rm $o/calls.unsorted.vcf
fi
for regions in clonal clonal_and_subclonal; do
  d=$O/bench_nocap/$regions; rm -rf $d; mkdir -p $d/vcf
  cp $R/bench/$regions/gt.vcf $d/gt.vcf
  $PYTHON $HELPERS/prep_vcf.py $o/calls.vcf $d/vcf/VESPA_t05.vcf --bed $V05/GRCh38_HG008-T-V0.5_somatic-stvar-$regions.draftbenchmark.bed --pass-only
  echo -e "$d/vcf/VESPA_t05.vcf\tVESPA_t05" > $d/callers.tsv
  $MINDA_PYTHON $MINDA truthset --base $d/gt.vcf --tsv $d/callers.tsv --out_dir $d/out --tolerance 500 --min_size 50 \
    --sample_name HG008 > $d/minda.log 2>&1
  $PYTHON $HELPERS/minda_excluded_truth.py $d $R/bench/truth_all.vcf $MINDA_PYTHON $MINDA
  { echo -e "method\tcallset\tTP\tFP\tFN"
    awk -F'\t' 'NR > 1 {print "minda\t" $1 "\t" $2 "\t" $4 "\t" $3}' $d/out/results/overall_rule.tsv; } > $d/summary.tsv
done
echo "done"
