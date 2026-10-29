#!/usr/bin/env bash
# Germline-removal cascade (Results, block A): rows 1, 1b and 2 of the table; rows 3-5 come from bp_coverage/
# and are benchmarked by bench.sh (already run on 2026-10-07).
# Order: 00 reads, 01 unitigs (hifiasm, ~4 h), PILOT search of chr20 against the reference-only index (sizes
# the SFS set before committing to the genome-wide run), 02 indexes and searches (~4-5 h), 03 calls and tables.
# Usage: nohup bash run_all.sh > logs/run_all.log 2>&1 &      (RUN_PILOT=1 runs the chr20 pilot; STOP_AFTER_PILOT=1 stops after it)
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"
HERE=$(cd "$(dirname "$0")" && pwd)
D=$HG008/reads_ec_align/germline_cascade
mkdir -p $D/logs
log() { echo "[$(date '+%F %T')] $*"; }
log 00_reads;   bash $HERE/00_reads.sh
log 01_unitigs; bash $HERE/01_unitigs.sh
grep -q IDENTICAL $D/unitigs/check_ec.txt || log "WARNING: regenerated normal EC reads differ from the paper run (see unitigs/check_ec.txt)"
RUN_PILOT=${RUN_PILOT:-0}   # the pilot only sized the SFS set; the calls are done for every row anyway
if [[ $RUN_PILOT == 1 ]]; then
log pilot_chr20; PILOT=chr20 INDEXES=ref bash $HERE/02_index_search.sh
# same reads against the full index of the paper: SFS of chr20 reads in search/tumor.dedup.cap500.sfs
P=$D/search_pilot_chr20
zcat $P/tumor.chr20.ec.fa.gz | awk '/^>/{print substr($1,2)}' | sort -u > $P/reads.txt
full=$(awk 'NR==FNR{r[$1]=1; next} /^#/{next} $1!="*"{cur=$1} (cur in r)' $P/reads.txt $HG008/reads_ec_align/search/tumor.dedup.cap500.sfs | wc -l)
log "pilot chr20: $(grep -vc '^#' $P/tumor.ref.dedup.cap500.sfs) SFS with the reference-only index, $full with the full index"
# pilot calls on the chr20 SFS, without and with the normal: time and memory to size the genome-wide calls
R=$HG008/reads_ec_align
for nb in - $R/align/normal.ec.hp.bam; do
  tag=nogl; norm=(); [[ $nb != - ]] && { tag=with_normal; norm=(--normal-contigs-bam $nb); }
  [[ -s $P/call_$tag.vcf ]] && continue
  /usr/bin/time -v -o $P/time_call_$tag.txt $VESPA call \
    --reference $R/ref/GRCh38_GIABv3.chr1-22X.fa --bam $R/align/tumor.ec.hp.bam --sfs $P/tumor.ref.dedup.cap500.sfs "${norm[@]}" \
    --bed-exclusion $REFS/grch38_centromers.bed --call-provenance $P/prov_$tag.tsv \
    --min-cluster-weight 5 --min-sv-length 50 --min-mapq 60 --clipped --no-clipped-fallback --threads 24 \
    > $P/call_$tag.vcf 2> $P/call_$tag.log
  log "pilot call $tag: $(grep -E 'Elapsed|Maximum resident' $P/time_call_$tag.txt | tr -s ' ' | tr '\n' ' ')"
done
fi
STOP_AFTER_PILOT=${STOP_AFTER_PILOT:-0}
[[ ${STOP_AFTER_PILOT:-0} == 1 ]] && { log "stopped after the pilot"; exit 0; }
log 02_index_search; bash $HERE/02_index_search.sh
log 03_call_tables;  R1_CALLS=${R1_CALLS:-1} bash $HERE/03_call_tables.sh
log done
