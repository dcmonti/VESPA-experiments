#!/usr/bin/env bash
# One cell: EC of the tumor of one sample with one setting, alignment + haplotag, search,
# call without the normal with the provenance dump + projection of the SFS, call with the normal at
# weights WEIGHTS.   USAGE: bash 05_cell.sh t90_p0.6 fnr0.10|classic|raw
# raw (2026-09-30): tumor without EC, control; raw reads and raw BAM already aligned (raw/, align_tag
# as for the other cells); the normal and the index remain the corrected ones, as in all cells.
set -euo pipefail
source "$(dirname "$0")/config.sh"; source $SRC/lib_sfs.sh
S=${1:?sample}; EC=${2:?setting}; RHO=${S#*_p}
W=$D/cells/$S/$EC; mkdir -p $W; N=$D/normal
case $EC in
  classic|raw) FLAG="" ;;
  fnr*)    FLAG="--tumor-purity $RHO --clonality $CLONALITY --somatic-fnr ${EC#fnr}" ;;
  *) echo "unknown setting: $EC" >&2; exit 1 ;;
esac
if [[ $EC == raw && ! -s $W/tumor.ec.fa.gz ]]; then
  log "$S raw: raw reads, no EC"
  zcat $D/reads/$S.fq.gz | awk 'NR%4==1{print ">"substr($1,2)} NR%4==2{print}' | gzip -1 > $W/tumor.ec.fa.gz
  ln -sf $D/raw/$S.bam $W/tumor.ec.bam; ln -sf $D/raw/$S.bam.bai $W/tumor.ec.bam.bai
  echo "no EC (raw reads)" > $W/threshold.txt
fi
if [[ ! -s $W/tumor.ec.fa.gz ]]; then
  log "$S $EC: hifiasm $FLAG"; rm -rf $W/work; mkdir -p $W/work
  /usr/bin/time -v $HIFIASM -t $THREADS --write-ec --bin-only $FLAG -o $W/work/ec $D/reads/$S.fq.gz > $W/hifiasm.log 2>&1
  gzip -1 -c $W/work/ec.ec.fa > $W/tumor.ec.fa.gz && rm -rf $W/work
  grep -a -E "NB somatic threshold|using user-specified het_cov|peak_hom" $W/hifiasm.log | tail -3 > $W/threshold.txt || true
fi
[[ -s $W/tumor.ec.bam.bai ]] || { log "$S $EC: alignment"; align_tag $W/tumor.ec.fa.gz $W/tumor.ec.bam; }
[[ -s $W/stats.txt ]] || samtools stats $W/tumor.ec.bam > $W/stats.txt
[[ -s $W/tumor.dedup.cap500.sfs ]] || { log "$S $EC: search"; search_sfs $N/index/full.fmd $W/tumor.ec.fa.gz $W/tumor; }
if [[ ! -s $W/prov/provenance.tsv ]]; then
  mkdir -p $W/prov; log "$S $EC: call without the normal with provenance"
  $VESPA call --reference $D/ref/region.fa --bam $W/tumor.ec.bam --sfs $W/tumor.dedup.cap500.sfs \
    --clip-clusters $W/prov/calls.clipped_clusters --call-provenance $W/prov/provenance.tsv \
    --threads $THREADS --min-cluster-weight 5 --min-sv-length 50 --min-mapq 60 --max-cluster-dist 0 \
    --clipped --no-clipped-fallback > $W/prov/calls.vcf 2> $W/prov/call.log
  python $S2G $W/tumor.dedup.cap500.sfs $W/tumor.ec.bam $W/prov/sfs_supp.bed.gz $THREADS 2> $W/prov/sfs_to_genome.log
  zcat $W/prov/sfs_supp.bed.gz | sort -k4,4 -k5,5n -k6,6n -S2G | awk -F'\t' '
    function flush(){ if(n==0) return; for(i=1;i<=n;i++) if(!hasexact || ex[i]==1) print buf[i] }
    { key=$4"\t"$5"\t"$6; if(key!=prev){ flush(); n=0; hasexact=0; prev=key }
      n++; buf[n]=$0; ex[n]=$7; if($7==1) hasexact=1 } END{ flush() }' | sort -k1,1 -k2,2n -S2G | gzip > $W/prov/sfs_clean.bed.gz
fi
if [[ ! -s $W/prov_w1/provenance.tsv ]]; then
  # m_min 1 only for the dump: in pcall a cluster below m_min before splitting writes no
  # EVAL rows, so at m_min 5 the "clusters of any weight" row would contain only clusters already >= 5
  mkdir -p $W/prov_w1; log "$S $EC: provenance dump at m_min 1"
  $VESPA call --reference $D/ref/region.fa --bam $W/tumor.ec.bam --sfs $W/tumor.dedup.cap500.sfs \
    --clip-clusters $W/prov_w1/calls.clipped_clusters --call-provenance $W/prov_w1/provenance.tsv \
    --threads $THREADS --min-cluster-weight 1 --min-sv-length 50 --min-mapq 60 --max-cluster-dist 0 \
    --clipped --no-clipped-fallback > $W/prov_w1/calls.vcf 2> $W/prov_w1/call.log
fi
for w in $WEIGHTS; do O=$W/call_w$w; [[ -s $O/calls.vcf ]] && continue; mkdir -p $O
  $VESPA call --reference $D/ref/region.fa --bam $W/tumor.ec.bam --sfs $W/tumor.dedup.cap500.sfs \
    --normal-contigs-bam $N/normal.ec.bam --clip-clusters $O/calls.clipped_clusters \
    --threads $THREADS --min-cluster-weight $w --min-sv-length 50 --min-mapq 60 --max-cluster-dist 0 \
    --clipped --no-clipped-fallback > $O/calls.unsorted.vcf 2> $O/call.log
  bcftools sort $O/calls.unsorted.vcf -o $O/calls.vcf 2>/dev/null && rm -f $O/calls.unsorted.vcf
done
log "$S $EC: done ($(cat $W/threshold.txt | tail -1))"
