#!/usr/bin/env bash
# HiFi reads of one sample: for each component sequence (hap1, hap2, hap1_T1, hap1_T2) PBSIM3
# simulates the subreads (QSHMM model, PASSES passes, lengths READ_LEN_MEAN +- READ_LEN_SD) at the
# depth of the component x YIELD, ccs builds their consensus (rq >= 0.99, as for real HiFi), and the
# reads are renamed <sample>.<component>.<n>: the name carries the haplotype of origin.
# Depth per component:
#   normal      hap1 NORMAL_DEPTH/2, hap2 NORMAL_DEPTH/2
#   tC_pRHO     hap2 C/2, hap1 C/2 (1-RHO), hap1_T1 C/2 RHO/2, hap1_T2 C/2 RHO/2
# Different seed per sample and component: the reads of the normal and the normal components of the tumor
# are independent.
# USAGE: bash 03_simulate.sh normal | t90_p0.6 | ...
set -euo pipefail
source "$(dirname "$0")/config.sh"
S=${1:?sample}; O=$D/reads/$S; mkdir -p $O; FQ=$D/reads/$S.fq.gz
[[ -s $FQ ]] && { log "$S: reads already simulated"; exit 0; }
QSHMM=$(find $SIMENV -name "QSHMM-RSII.model" | head -1); [[ -s $QSHMM ]] || { echo "QSHMM model missing (00_env.sh)" >&2; exit 1; }
if [[ $S == normal ]]; then
  COMPS="hap1:$(awk -v c=$NORMAL_DEPTH 'BEGIN{print c/2}') hap2:$(awk -v c=$NORMAL_DEPTH 'BEGIN{print c/2}')"
else
  C=${S#t}; C=${C%%_*}; RHO=${S#*_p}
  COMPS=$(awk -v c=$C -v r=$RHO 'BEGIN{printf "hap2:%g hap1:%g hap1_T1:%g hap1_T2:%g", c/2, c/2*(1-r), c/2*r/2, c/2*r/2}')
fi
log "$S: components $COMPS"
i=0; : > $O/depths.tsv
for cd in $COMPS; do
  comp=${cd%%:*}; dep=${cd#*:}; i=$((i+1))
  awk -v d=$dep 'BEGIN{exit !(d>0)}' || continue
  [[ -s $O/$comp.fq.gz ]] && continue
  seed=$(( SEED + $(echo "$S.$comp" | cksum | cut -d' ' -f1) % 100000 ))
  log "$S $comp: PBSIM3 at ${dep}x (x $YIELD), seed $seed"
  ( cd $O && $SIMENV/bin/pbsim --strategy wgs --method qshmm --qshmm $QSHMM \
      --depth $(awk -v d=$dep -v y=$YIELD 'BEGIN{print d*y}') --genome $D/ref/$comp.fa \
      --length-mean $READ_LEN_MEAN --length-sd $READ_LEN_SD --pass-num $PASSES --seed $seed \
      --prefix $comp > $comp.pbsim.log 2>&1 )
  # PBSIM3 >= 3.0.5 writes the subreads as BAM (<prefix>_0001.bam); earlier versions as SAM
  if [[ -s $O/${comp}_0001.bam ]]; then mv $O/${comp}_0001.bam $O/$comp.subreads.bam
  elif [[ -s $O/${comp}_0001.sam ]]; then $SIMENV/bin/samtools view -@ 4 -b -o $O/$comp.subreads.bam $O/${comp}_0001.sam
  else echo "$S $comp: PBSIM3 did not write the subreads (see $O/$comp.pbsim.log)" >&2; exit 1; fi
  rm -f $O/${comp}_0001.sam $O/${comp}_0001.maf* $O/${comp}_0001.ref
  log "$S $comp: ccs"
  $SIMENV/bin/ccs $O/$comp.subreads.bam $O/$comp.ccs.fq.gz -j $THREADS --min-rq 0.99 > $O/$comp.ccs.log 2>&1
  rm -f $O/$comp.subreads.bam*
  zcat $O/$comp.ccs.fq.gz | awk -v p="$S.$comp" 'NR%4==1{n++; print "@"p"."n; next} {print}' | gzip -1 > $O/$comp.fq.gz
  rm -f $O/$comp.ccs.fq.gz*
  L=$(samtools faidx $D/ref/$comp.fa 2>/dev/null; cut -f2 $D/ref/$comp.fa.fai)
  zcat $O/$comp.fq.gz | awk -v c=$comp -v t=$dep -v L=$L 'NR%4==2{n++; b+=length($0)} END{printf "%s\t%g\t%d\t%.2f\n", c, t, n, b/L}' >> $O/depths.tsv
done
cat $(for cd in $COMPS; do c=${cd%%:*}; [[ -s $O/$c.fq.gz ]] && echo $O/$c.fq.gz; done) > $FQ
log "$S: obtained depths (component, target, reads, x):"; cat $O/depths.tsv
