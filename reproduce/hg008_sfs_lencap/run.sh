#!/usr/bin/env bash
# HG008: effect of the SFS length cap (dedup_sfs.py --max-length 500 in 03_vespa_index_search.sh).
# Same run as reads_ec_align (index, reads, alignments, calls, bench); the only change: the tumor SFS are not capped.
#   1. tumor search on the full germline index, as 03_vespa_index_search.sh (the corrected tumor reads are those
#      recovered from the EC BAM by germline_cascade/00_reads.sh)
#   2. dedup_sfs.py without a cap (tumor.dedup.nocap.sfs) and with the 500 bp cap (tumor.dedup.cap500.sfs, a check:
#      it must equal reads_ec_align/search/tumor.dedup.cap500.sfs)
#   3. VESPA call on the uncapped SFS, only s = 5 (the default run of 04_call.sh), 32 threads (as the default run
#      of 04_call.sh, so that time and memory compare with reads_ec_align/time/vespa.t05.txt)
#   4. Minda bench of the VESPA runs, as 05_bench.sh (PASS, regions, excluded-truth rule), clonal and c+s
# Output: $O/{search,calls,bench,time,logs}. Usage: bash run.sh
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"
R=$HG008/reads_ec_align
O=$R/sfs_lencap
HERE=$(cd "$(dirname "$0")" && pwd)
HELPERS=$(cd $HERE/../bench && pwd)
V05=$HG008/ground_truths/v0.5
REF=$R/ref/GRCh38_GIABv3.chr1-22X.fa
EXCLUDE=$REFS/grch38_centromers.bed
TUMOR_EC=$R/germline_cascade/reads/tumor.ec.fa.gz
THRESHOLDS="5"   # the weight of the analyses of S3 (2026-10-15); the grid is not needed here
mkdir -p $O/search $O/calls $O/bench $O/time $O/logs
timed() { local name=$1; shift; /usr/bin/time -v -o $O/time/$name.txt "$@"; }
export PATH=$SV_BIN:$PATH

# 1. search
if [[ ! -s $O/search/tumor.sfs ]]; then
  timed vespa_search_tumor $VESPA search --index $R/index/full.fmd --fastx $TUMOR_EC --bsize 10000 --threads 48 \
    --hpc --hpc-cap 2 > $O/search/tumor.sfs.part 2> $O/logs/vespa_search_tumor.log
  grep -qa "All done!" $O/logs/vespa_search_tumor.log || { echo "search failed" >&2; exit 1; }
  mv $O/search/tumor.sfs.part $O/search/tumor.sfs
fi
# 2. dedup, without and with the cap
[[ -s $O/search/tumor.dedup.nocap.sfs ]] || timed dedup_nocap $PYTHON $LIB/dedup_sfs.py $O/search/tumor.sfs \
  $O/search/tumor.dedup.nocap.sfs --max-length 0 2> $O/logs/dedup_nocap.log
[[ -s $O/search/tumor.dedup.cap500.sfs ]] || timed dedup_cap500 $PYTHON $LIB/dedup_sfs.py $O/search/tumor.sfs \
  $O/search/tumor.dedup.cap500.sfs --max-length 500 2> $O/logs/dedup_cap500.log
if cmp -s $O/search/tumor.dedup.cap500.sfs $R/search/tumor.dedup.cap500.sfs; then echo "cap500 check: identical" > $O/logs/cap500_check.txt
else echo "cap500 check: DIFFERENT from reads_ec_align/search/tumor.dedup.cap500.sfs" > $O/logs/cap500_check.txt; fi
# 3. calls
for t in $THRESHOLDS; do
  run=$(printf "t%02d" $t); o=$O/calls/$run
  [[ -s $o/calls.vcf ]] && continue
  rm -rf $o; mkdir -p $o
  timed vespa.$run $VESPA call --reference $REF --bam $R/align/tumor.ec.hp.bam --sfs $O/search/tumor.dedup.nocap.sfs \
    --normal-contigs-bam $R/align/normal.ec.hp.bam --bed-exclusion $EXCLUDE --min-cluster-weight $t --min-sv-length 50 \
    --min-mapq 60 --clipped --no-clipped-fallback --threads 32 > $o/calls.unsorted.vcf 2> $o/call.log
  $BCFTOOLS sort $o/calls.unsorted.vcf -o $o/calls.vcf && rm $o/calls.unsorted.vcf
done
# 4. bench (Minda, PASS, regions, excluded-truth rule)
TRUTH_ALL=$R/bench/truth_all.vcf
[[ -s $TRUTH_ALL ]] || $PYTHON $HELPERS/prep_vcf.py $V05/GRCh38_HG008-T-V0.5_somatic-stvar_ALL.draftbenchmark.vcf $O/bench/truth_all.vcf --min-len 50
[[ -s $TRUTH_ALL ]] || TRUTH_ALL=$O/bench/truth_all.vcf
for regions in clonal clonal_and_subclonal; do
  d=$O/bench/$regions; rm -rf $d; mkdir -p $d/vcf
  cp $R/bench/$regions/gt.vcf $d/gt.vcf
  : > $d/callers.tsv
  for t in $THRESHOLDS; do
    name=$(printf "VESPA_t%02d" $t)
    $PYTHON $HELPERS/prep_vcf.py $O/calls/$(printf "t%02d" $t)/calls.vcf $d/vcf/$name.vcf --bed $V05/GRCh38_HG008-T-V0.5_somatic-stvar-$regions.draftbenchmark.bed --pass-only
    echo -e "$d/vcf/$name.vcf\t$name" >> $d/callers.tsv
  done
  $MINDA_PYTHON $MINDA truthset --base $d/gt.vcf --tsv $d/callers.tsv --out_dir $d/out --tolerance 500 --min_size 50 \
    --sample_name HG008 > $d/minda.log 2>&1
  $PYTHON $HELPERS/minda_excluded_truth.py $d $TRUTH_ALL $MINDA_PYTHON $MINDA
  { echo -e "method\tcallset\tTP\tFP\tFN"
    awk -F'\t' 'NR > 1 {print "minda\t" $1 "\t" $2 "\t" $4 "\t" $3}' $d/out/results/overall_rule.tsv; } > $d/summary.tsv
done
echo "done"
