#!/usr/bin/env bash
# Breakpoint coverage by stage (tab:res-bp-coverage, Results 3.3) on the reads_ec_align HG008 run.
# Same procedure as het6/hpc_caps (bp_point_table_prov.py, TOL 50), on the current inputs:
#   SFS    reads_ec_align/search/tumor.dedup.cap500.sfs
#   BAM    reads_ec_align/align/tumor.ec.hp.bam (EC reads, minimap2 map-hifi default, GIABv3 chr1-22X)
#   VESPA  tools/vespa/VESPA (md5 669e3bed, has --call-provenance), flags of 04_call.sh, w = 5
# Steps:
#   1. project the SFS on the genome (sfs_to_genome.py, supplementary alignments included), keep the
#      exact placements when a SFS has one -> sfs_clean.bed.gz          (placed SFS row)
#   2. VESPA call WITHOUT the normal, --call-provenance                  (clusters, candidate clusters)
#   3. VESPA call WITH the normal, --call-provenance                     (check only: every breakpoint
#      covered by a PASS call must be covered by a candidate cluster)
#   4. table, clonal and clonal_and_subclonal regions of GIAB V0.5, tolerance 50 bp
# Usage: bash run.sh
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"
export PATH=$SV_BIN:$PATH
R=$HG008/reads_ec_align
O=$R/bp_coverage
SFS=$R/search/tumor.dedup.cap500.sfs
BAM=$R/align/tumor.ec.hp.bam
NBAM=$R/align/normal.ec.hp.bam
REF=$R/ref/GRCh38_GIABv3.chr1-22X.fa
EXCL=$REFS/grch38_centromers.bed
T=${THREADS:-24}

# 1. placed SFS
if [[ ! -s $O/sfs_clean.bed.gz ]]; then
  python $LIB/sfs_to_genome.py $SFS $BAM $O/sfs_supp.bed.gz $T 2> $O/sfs_to_genome.log
  zcat $O/sfs_supp.bed.gz | sort -k4,4 -k5,5n -k6,6n -S4G | awk -F'\t' '
    function flush(){ if(n==0) return; for(i=1;i<=n;i++) if(!hasexact || ex[i]==1) print buf[i] }
    { key=$4"\t"$5"\t"$6
      if(key!=prev){ flush(); n=0; hasexact=0; prev=key }
      n++; buf[n]=$0; ex[n]=$7; if($7==1) hasexact=1 }
    END{ flush() }' | sort -k1,1 -k2,2n -S4G | gzip > $O/sfs_clean.bed.gz.part
  mv $O/sfs_clean.bed.gz.part $O/sfs_clean.bed.gz
fi

# 2-3. calls with provenance, without and with the normal
call() { # outdir normal_bam|-
  local d=$O/$1 nb=$2; mkdir -p $d
  [[ -s $d/calls.vcf ]] && return 0
  local norm=(); [[ $nb != - ]] && norm=(--normal-contigs-bam $nb)
  $VESPA call --reference $REF --bam $BAM --sfs $SFS "${norm[@]}" --bed-exclusion $EXCL \
    --clip-clusters $d/calls.clipped_clusters --call-provenance $d/provenance.tsv \
    --min-cluster-weight 5 --min-sv-length 50 --min-mapq 60 --clipped --no-clipped-fallback --threads $T \
    > $d/calls.unsorted.vcf 2> $d/call.log
  bcftools sort $d/calls.unsorted.vcf -o $d/calls.vcf 2>/dev/null && rm -f $d/calls.unsorted.vcf
}
call nogl_w5 - &
call with_normal_w5 $NBAM &
wait

# 4. table (and the same table on the run with the normal, for the check)
FAI=$REF.fai TOL=50 python $REPRO/bp_coverage/bp_point_table_prov.py $O $O/nogl_w5/provenance.tsv HG008_nogl_w5 > $O/bp_point_prov_tol50.tsv
FAI=$REF.fai TOL=50 python $REPRO/bp_coverage/bp_point_table_prov.py $O $O/with_normal_w5/provenance.tsv HG008_withnormal_w5 > $O/bp_point_prov_withnormal_tol50.tsv
echo done
