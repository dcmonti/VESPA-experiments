#!/usr/bin/env bash
# HG002, chr1: assign each read to the homolog it comes from with the TRUTH phasing, not with ours, so that q is
# measured apart from the errors of the haplotype tags (assumption (iii) of the model).
#
#   phased heterozygous SNVs of the dipcall VCF, inside the small-variant benchmark regions, all in one phase set
#   (the two assembly haplotypes are chromosome-scale) -> whatshap haplotag. HP=1 is the first allele of GT (hap1 of
#   the assembly), HP=2 the second; a read with no informative SNV stays untagged.
#
# Output, in $OUT: truth/hg002.<chrom>.snv.phased.vcf.gz, align/hg002.$TAG.raw.hp.bam(.bai), logs/haplotag*.log
# Usage: bash 02_haplotag.sh            (chr1)
#        CHROMS="chr1 ... chr22" TAG=autosomes bash 02_haplotag.sh   (one whatshap per chromosome, merged)
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"

OUT=${OUT:-$HG002/q_visibility}
REF=$HG008/reads_ec_align/ref/GRCh38_GIABv3.chr1-22X.fa
DIP=$OUT/truth/GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.vcf.gz
SMVAR_BED=$OUT/truth/HG002_GRCh38_v5.0q_smvar.benchmark.bed
CHROMS=${CHROMS:-chr1}
TAG=${TAG:-chr1}
JOBS=${JOBS:-11}

phase_vcf() { # chrom -> truth/hg002.<chrom>.snv.phased.vcf.gz
  local c=$1 vcf=$OUT/truth/hg002.$1.snv.phased.vcf.gz
  [[ -s $vcf.tbi ]] && return 0
  # PS added explicitly (one block per chromosome), so that the phase set does not depend on how whatshap reads a
  # phased GT without PS
  $SV_BIN/bcftools view -r $c -T $SMVAR_BED -v snps -g het -p $DIP \
    | awk 'BEGIN{OFS="\t"}
           /^##/ {print; next}
           /^#CHROM/ {print "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set: one per chromosome, the haplotypes of the diploid assembly\">"; print; next}
           {$9 = $9 ":PS"; $10 = $10 ":1"; print}' \
    | $SV_BIN/bgzip > $vcf.part
  mv $vcf.part $vcf
  $SV_BIN/tabix -p vcf $vcf
}

tag_chrom() { # chrom -> align/hg002.$TAG.<chrom>.hp.bam (for TAG=chr1: the final file directly)
  local c=$1 bam=$OUT/align/hg002.$TAG.raw.bam vcf=$OUT/truth/hg002.$1.snv.phased.vcf.gz
  local hp=$OUT/align/hg002.$TAG.$c.hp.bam log=$OUT/logs/haplotag.$TAG.$c.log
  [[ $TAG == chr1 ]] && { hp=$OUT/align/hg002.chr1.raw.hp.bam; log=$OUT/logs/haplotag.log; }
  [[ -s $hp.bai || -s $OUT/align/hg002.$TAG.raw.hp.bam.bai ]] && return 0
  local sample=$($SV_BIN/bcftools query -l $vcf)
  $SV_BIN/whatshap haplotag --reference $REF --ignore-read-groups --sample $sample --regions $c \
    --output-haplotag-list $OUT/align/hg002.$TAG.$c.haplotag.tsv.gz -o $hp.part $vcf $bam 2> $log
  mv $hp.part $hp
  $SV_BIN/samtools index $hp
}
export -f phase_vcf tag_chrom
export BIN OUT REF DIP SMVAR_BED TAG

for c in $CHROMS; do phase_vcf $c; done
printf "%s\n" $CHROMS | xargs -P $JOBS -I{} bash -c 'tag_chrom {}'
hp=$OUT/align/hg002.$TAG.raw.hp.bam
if [[ $TAG != chr1 && ! -s $hp.bai ]]; then
  $SV_BIN/samtools merge -@ 16 -f -o $hp.part $(for c in $CHROMS; do echo $OUT/align/hg002.$TAG.$c.hp.bam; done)
  mv $hp.part $hp
  $SV_BIN/samtools index -@ 16 $hp
  rm -f $OUT/align/hg002.$TAG.chr*.hp.bam $OUT/align/hg002.$TAG.chr*.hp.bam.bai
fi
echo done
