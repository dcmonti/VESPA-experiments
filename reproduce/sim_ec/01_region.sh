#!/usr/bin/env bash
# Reference = hg38 region (sequence "$CTG"); hap1 and hap2 = the corresponding segments of the verkko,
# found by aligning the hg38 region to the diploid assembly (minimap2 -x asm5). Then hap1 and hap2
# are realigned to the reference (--cs=long): the hap1 -> reference map is used to write the truth
# in reference coordinates, and the germline differences to keep the somatic variants away from them.
set -euo pipefail
source "$(dirname "$0")/config.sh"; mkdir -p $D/ref $D/logs
R=$D/ref/region.fa
if [[ ! -s $R.fai ]]; then
  samtools faidx $HG38 $REGION | sed "1s/.*/>$CTG/" > $R && samtools faidx $R
  log "region $REGION: $(cut -f2 $R.fai) bp, N: $(grep -v '>' $R | tr -cd 'Nn' | wc -c)"
fi
P=$D/ref/region_vs_verkko.paf
# secondary alignments too (-N 10): relative to hap1, the alignment to hap2 is a secondary
[[ -s $P ]] || minimap2 -x asm5 -c -t $THREADS --secondary=yes -N 10 $VERKKO $R > $P 2> $D/logs/region_vs_verkko.log
for h in hap1 hap2; do
  [[ -s $D/ref/$h.fa ]] && continue
  # contig and strand with the most aligned bases; extent = min/max of the coordinates on the contig
  read -r ctg strand ts te qcov < <(awk -v h=$h -F'\t' '$6 ~ "^"h"_" && $12 >= 0 {k=$6"\t"$5; m[k]+=$10; if(!(k in s)||$8<s[k])s[k]=$8; if($9>e[k])e[k]=$9; qc[k]+=$4-$3; ql=$2}
      END{for(k in m) if(m[k]>best){best=m[k]; bk=k}; split(bk,a,"\t"); printf "%s %s %d %d %.3f\n", a[1], a[2], s[bk], e[bk], (ql>0?qc[bk]/ql:0)}' $P)
  [[ -n ${ctg:-} && ${ctg:-} != "" ]] || { echo "$h: no alignment of the region" >&2; exit 1; }
  log "$h: $ctg $strand $ts-$te, fraction of the region aligned $qcov"
  awk -v q=$qcov 'BEGIN{exit !(q>=0.95)}' || { echo "$h: less than 95% of the region aligned to a single contig, choose another REGION" >&2; exit 1; }
  samtools faidx $VERKKO "$ctg:$((ts+1))-$te" > $D/ref/$h.tmp.fa
  if [[ $strand == "-" ]]; then samtools faidx -i $D/ref/$h.tmp.fa "$ctg:$((ts+1))-$te" > $D/ref/$h.rc.fa; mv $D/ref/$h.rc.fa $D/ref/$h.tmp.fa; fi
  sed "1s/.*/>$h/" $D/ref/$h.tmp.fa > $D/ref/$h.fa && rm $D/ref/$h.tmp.fa && samtools faidx $D/ref/$h.fa
done
for h in hap1 hap2; do
  [[ -s $D/ref/${h}_vs_ref.paf ]] || minimap2 -x asm5 -c --cs=long -t $THREADS --secondary=no $R $D/ref/$h.fa > $D/ref/${h}_vs_ref.paf 2>> $D/logs/hap_vs_ref.log
  log "$h -> reference: $(awk -F'\t' '$12>=0{n++; m+=$10; b+=$11} END{printf "%d alignments, identity %.4f", n, m/b}' $D/ref/${h}_vs_ref.paf)"
done
