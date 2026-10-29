#!/usr/bin/env bash
# Germline-removal cascade, step 1: the unitigs of the normal (hifiasm r_utg), regenerated.
# 03_vespa_index_search.sh deleted ec/normal.r_utg.fa.gz. They are rebuilt exactly as in
# reproduce/hg008/01_reads_ec_align.sh: normal reads with primary alignment on chr1-22, X of the NIST BAM,
# one hifiasm run per chromosome (v1.0.0-somatic, --het-cov 4, full run), unitigs prefixed by chromosome.
# Check: the corrected normal reads of this run are compared with those of the paper run (00_reads.sh);
# identical reads mean the hifiasm runs are reproduced, hence the same unitigs.
# Output: unitigs/normal.r_utg.fa.gz, unitigs/check_ec.txt
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"
R=$HG008/reads_ec_align
D=$R/germline_cascade
NORMAL_BAM=$HG008/jan2024_2026-09-20/raw/HG008-N-P.jan2024.bam
THREADS=48
CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"
U=$D/unitigs
mkdir -p $U $D/time $D/logs
[[ -s $U/normal.r_utg.fa.gz ]] && { echo "already done"; exit 0; }
for chrom in $CHROMS; do
  [[ -s $U/normal.$chrom.r_utg.fa.gz ]] && continue
  fq=$U/normal.$chrom.fq.gz; work=$U/normal.$chrom.work
  $SAMTOOLS view -u -F 0x900 $NORMAL_BAM $chrom | $SAMTOOLS fastq - 2>/dev/null | gzip -1 > $fq
  rm -rf $work; mkdir -p $work
  /usr/bin/time -v -o $D/time/hifiasm.normal.$chrom.txt $HIFIASM -t $THREADS --write-ec --het-cov 4 -o $work/ec $fq \
    > $D/logs/hifiasm.normal.$chrom.log 2>&1
  awk -v p=${chrom}_ '/^S/ {print ">" p $2; print $3}' $work/ec.bp.r_utg.gfa | gzip -1 > $U/normal.$chrom.r_utg.fa.gz.part
  gzip -1 -c $work/ec.ec.fa > $U/normal.$chrom.ec.fa.gz
  mv $U/normal.$chrom.r_utg.fa.gz.part $U/normal.$chrom.r_utg.fa.gz
  rm -rf $work $fq
done
cat $(for c in $CHROMS; do echo $U/normal.$c.r_utg.fa.gz; done) > $U/normal.r_utg.fa.gz.part
mv $U/normal.r_utg.fa.gz.part $U/normal.r_utg.fa.gz

# check: same corrected reads as the paper run (name and sequence, order-independent)
digest() { zcat "$@" | awk '/^>/{n=substr($1,2); next} {print n"\t"$0}' | sort -S4G | md5sum | cut -d' ' -f1; }
a=$(digest $(for c in $CHROMS; do echo $U/normal.$c.ec.fa.gz; done)); b=$(digest $D/reads/normal.ec.fa.gz)
{ echo "regenerated $a"; echo "paper_run   $b"; [[ $a == $b ]] && echo "IDENTICAL" || echo "DIFFERENT"; } > $U/check_ec.txt
cat $U/check_ec.txt
rm -f $U/normal.chr*.ec.fa.gz $U/normal.chr*.r_utg.fa.gz
