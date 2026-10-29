#!/usr/bin/env bash
# Germline-removal cascade, step 3: calls, breakpoint coverage and Minda for the new rows.
# Table layout (2026-10-08): one row per index, calls in three columns: no normal, in-call filter (all calls of the
# run with the normal), in-call and post-call filters (PASS calls of the same run).
#   R1_ref_nofilters     reference-only index, VESPA call without the normal          -> stages + calls "no normal"
#   R1b_ref_normal       reference-only index, VESPA call with the normal             -> R1b_ref_incall, R1b_ref_pass
#   R2_refutg_nofilters  reference + unitigs index, VESPA call without the normal     -> stages + calls "no normal"
#   R2b_refutg_normal    reference + unitigs index, VESPA call with the normal        -> R2b_refutg_incall, R2b_refutg_pass
# The normal acts only on the candidate SVs, so the runs with the normal get no breakpoint table (their stages are
# those of the run without it, same index and SFS).
# Calls: VESPA with the flags of reproduce/hg008/04_call.sh, w = 5, --call-provenance (as bp_coverage/).
# Breakpoint coverage: SFS projected on the EC tumor BAM, then bp_point_table_prov.py (tolerance 50 bp;
# n, area and candidate SVs inside the regions in the last columns).
# Minda: the rows are appended to rows.tsv and bench.sh is run again (all rows, excluded-truth rule).
# Output: rows/<row>/{calls.vcf,provenance.tsv,sfs_clean.bed.gz,bp_point_prov_tol50.tsv}, bench/summary.tsv
#
# Node safety (2026-10-08; the reference-only index gives 225.9 M SFS, ~190x the full germline set, while other
# pipelines share the node):
#   - order: the reference + unitigs rows first, with their Minda, then the reference-only rows, so that a failure
#     or a very long run on the reference-only index does not hold back the other rows;
#   - every step of a row may fail (memory cap, timeout) without stopping the script: the row is reported as failed
#     and the next one goes on;
#   - memory: address-space cap per process (ulimit -v, MEM_CAP_GB). The VESPA call is one process, so the cap holds
#     for the whole call, including the hierarchical agglomerative clustering (n x n matrix and merge trace per
#     interval). sfs_to_genome.py forks its workers after loading every SFS in memory, and copy-on-write pages grow
#     in each worker: the cap does not bound their sum, so the reference-only projection runs with PROJ_THREADS_R1
#     workers (default 4) instead of THREADS, as is any other SFS file over 1 GB;
#   - time: each VESPA call runs under timeout CALL_TIMEOUT (default 12h; the HAC is O(n^3) per interval, and the
#     reference-only intervals hold many more SFS).
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"
export PATH=$SV_BIN:$PATH
R=$HG008/reads_ec_align
D=$R/germline_cascade
REF=$R/ref/GRCh38_GIABv3.chr1-22X.fa
BAM=$R/align/tumor.ec.hp.bam
NBAM=$R/align/normal.ec.hp.bam
EXCL=$REFS/grch38_centromers.bed
T=${THREADS:-24}
PROJ_THREADS_R1=${PROJ_THREADS_R1:-4}
CALL_TIMEOUT=${CALL_TIMEOUT:-12h}
MEM_CAP_GB=${MEM_CAP_GB:-100}
ulimit -v $((MEM_CAP_GB * 1024 * 1024))

warn() { echo "[$(date '+%F %T')] $*" >&2; }

# project SFS on the tumor EC BAM, keep the exact placements of each SFS when there is one: sfs workers out
project() {
  local sfs=$1 nt=$2 d=$3
  python $LIB/sfs_to_genome.py $sfs $BAM $d/sfs_supp.bed.gz $nt 2> $d/sfs_to_genome.log \
  && { zcat $d/sfs_supp.bed.gz | sort -k4,4 -k5,5n -k6,6n -S8G | awk -F'\t' '
      function flush(){ if(n==0) return; for(i=1;i<=n;i++) if(!hasexact || ex[i]==1) print buf[i] }
      { key=$4"\t"$5"\t"$6
        if(key!=prev){ flush(); n=0; hasexact=0; prev=key }
        n++; buf[n]=$0; ex[n]=$7; if($7==1) hasexact=1 }
      END{ flush() }' | sort -k1,1 -k2,2n -S8G | gzip > $d/sfs_clean.bed.gz.part; } \
  && mv $d/sfs_clean.bed.gz.part $d/sfs_clean.bed.gz \
  && rm -f $d/sfs_supp.bed.gz
}

row() { # name sfs normal_bam|- [call 1|0] [breakpoint table 1|0] [projection workers]
  local name=$1 sfs=$2 nb=$3 docall=${4:-1} dobp=${5:-1} nt=${6:-$T} d=$D/rows/$1; mkdir -p $d
  # any SFS file over 1 GB (the full germline set is 51 MB) is projected with the reduced number of workers
  [[ $(stat -c %s $sfs) -gt 1000000000 ]] && nt=$(( nt < PROJ_THREADS_R1 ? nt : PROJ_THREADS_R1 ))
  warn "$name: start"
  if [[ $docall == 0 ]]; then : > $d/provenance.tsv; fi   # SFS only: no clusters, no calls
  if [[ $docall == 1 && ! -s $d/calls.vcf ]]; then
    local norm=(); [[ $nb != - ]] && norm=(--normal-contigs-bam $nb)
    if ! /usr/bin/time -v -o $d/time_call.txt timeout $CALL_TIMEOUT $VESPA call --reference $REF --bam $BAM --sfs $sfs "${norm[@]}" \
         --bed-exclusion $EXCL --clip-clusters $d/calls.clipped_clusters --call-provenance $d/provenance.tsv \
         --min-cluster-weight 5 --min-sv-length 50 --min-mapq 60 --clipped --no-clipped-fallback --threads $T \
         > $d/calls.unsorted.vcf 2> $d/call.log; then
      warn "$name: VESPA call failed (memory cap ${MEM_CAP_GB} GB or timeout $CALL_TIMEOUT; exit status in $d/time_call.txt), see $d/call.log"
      # the SFS stage does not need the call: go on with an empty provenance (SFS columns only, as with call 0);
      # a later run of this script retries the call, since calls.vcf is missing
      mv $d/provenance.tsv $d/provenance.failed.tsv 2>/dev/null || true
      rm -f $d/calls.unsorted.vcf; : > $d/provenance.tsv; docall=0
    fi
    if [[ $docall == 1 ]]; then
      bcftools sort $d/calls.unsorted.vcf -o $d/calls.vcf 2>/dev/null && rm -f $d/calls.unsorted.vcf \
        || { warn "$name: bcftools sort failed"; return 0; }
    fi
  fi
  [[ $dobp == 1 ]] || { warn "$name: done (no breakpoint table)"; return 0; }
  if [[ ! -s $d/sfs_clean.bed.gz ]]; then
    project $sfs $nt $d || { warn "$name: SFS projection failed ($nt workers), see $d/sfs_to_genome.log"; rm -f $d/sfs_clean.bed.gz.part; return 0; }
  fi
  [[ $docall == 1 && ! -s $d/calls.vcf ]] && return 0   # call failed
  FAI=$REF.fai TOL=50 python $REPRO/bp_coverage/bp_point_table_prov.py $d $d/provenance.tsv $name > $d/bp_point_prov_tol50.tsv.part \
    && mv $d/bp_point_prov_tol50.tsv.part $d/bp_point_prov_tol50.tsv \
    || { warn "$name: breakpoint table failed"; rm -f $d/bp_point_prov_tol50.tsv.part; return 0; }
  warn "$name: done"
}
add_bench() { # bench_name run_name passonly
  [[ -s $D/rows/$2/calls.vcf ]] || return 0
  grep -q "^$1 " $D/rows.tsv || echo "$1    $D/rows/$2/calls.vcf    $3" >> $D/rows.tsv
}
bench() { bash $REPRO/germline_cascade/bench.sh || warn "bench.sh failed"; }

# R1_CALLS=0: the reference-only index gives too many SFS to call; row 1 keeps only the SFS stage, row 1b is dropped
R1_CALLS=${R1_CALLS:-1}
U=$D/search/tumor.ref_utg.dedup.cap500.sfs; F=$D/search/tumor.ref.dedup.cap500.sfs

# 1. reference + unitigs, then its Minda
row R2_refutg_nofilters $U -     1 1
row R2b_refutg_normal   $U $NBAM 1 0
add_bench R2_refutg_nofilters R2_refutg_nofilters 0
add_bench R2b_refutg_incall   R2b_refutg_normal   0
add_bench R2b_refutg_pass     R2b_refutg_normal   1
bench

# 2. reference only, then Minda again (bench.sh runs all rows of rows.tsv)
row R1_ref_nofilters    $F -     $R1_CALLS 1 $PROJ_THREADS_R1
[[ $R1_CALLS == 1 ]] && row R1b_ref_normal $F $NBAM 1 0
if [[ $R1_CALLS == 1 ]]; then
  add_bench R1_ref_nofilters R1_ref_nofilters 0
  add_bench R1b_ref_incall   R1b_ref_normal   0
  add_bench R1b_ref_pass     R1b_ref_normal   1
  bench
fi
for r in R1_ref_nofilters R2_refutg_nofilters; do [[ -s $D/rows/$r/bp_point_prov_tol50.tsv ]] && cat $D/rows/$r/bp_point_prov_tol50.tsv || true; done
