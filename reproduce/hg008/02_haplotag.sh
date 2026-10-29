#!/usr/bin/env bash
# HG008: haplotagging of the raw and error-corrected alignments written by 01_reads_ec_align.sh.
#
# Steps:
#   1. call small variants on the normal raw alignments with DeepVariant (PACBIO model)
#   2. phase them with whatshap phase, on the same alignments
#   3. assign each tumor and normal read to a haplotype with whatshap haplotag, on the raw alignments;
#      the assignments are saved as read lists (read name, haplotype, phase set, chromosome)
#   4. write HP and PS tags from the read lists into the raw and the error-corrected alignments, by read name
#      (error-corrected reads keep the name of the raw read they come from)
#
# Output, in $OUT:
#   haplotag/normal.phased.vcf.gz              DeepVariant calls, phased; the unphased DeepVariant VCF is deleted
#   haplotag/{tumor,normal}.haplotags.tsv      read lists
#   align/{tumor,normal}.{raw,ec}.hp.bam       tagged alignments, indexed; the untagged ones are deleted
#   time/<step>.txt                            GNU time (/usr/bin/time -v) report of each step
#
# Usage: bash 02_haplotag.sh
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"

OUT=${OUT:-$HG008/reads_ec_align}
REF=$OUT/ref/GRCh38_GIABv3.chr1-22X.fa

DEEPVARIANT=docker://google/deepvariant:1.9.0
SINGULARITY=/usr/bin/singularity
TAG_BAM=$LIB/tag_bam_from_list.py

THREADS=48

mkdir -p $OUT/haplotag $OUT/logs $OUT/time

# timed NAME COMMAND...: run COMMAND under GNU time, report in $OUT/time/NAME.txt.
timed() { local name=$1; shift; /usr/bin/time -v -o $OUT/time/$name.txt "$@"; }

DV_VCF=$OUT/haplotag/normal.dv.vcf.gz
PHASED=$OUT/haplotag/normal.phased.vcf.gz
if [[ ! -s $PHASED.tbi ]]; then
  # 1. DeepVariant on the normal.
  if [[ ! -s $DV_VCF.tbi ]]; then
    mkdir -p $OUT/haplotag/tmp
    timed deepvariant $SINGULARITY exec --cleanenv -B $OUT --env TMPDIR=$OUT/haplotag/tmp $DEEPVARIANT \
      run_deepvariant --model_type=PACBIO --ref=$REF --reads=$OUT/align/normal.raw.bam \
        --output_vcf=$DV_VCF --num_shards=$THREADS --intermediate_results_dir=$OUT/haplotag/dv_intermediate \
        > $OUT/logs/deepvariant.log 2>&1
    rm -rf $OUT/haplotag/tmp $OUT/haplotag/dv_intermediate
  fi

  # 2. Phasing on the normal. The phased VCF holds every DeepVariant record, so the unphased one is deleted.
  timed whatshap_phase $WHATSHAP phase --reference $REF --ignore-read-groups -o $PHASED $DV_VCF $OUT/align/normal.raw.bam \
    > $OUT/logs/whatshap_phase.log 2>&1
  timed tabix_phased $TABIX -f -p vcf $PHASED
  rm -f $OUT/haplotag/normal.dv.*
fi

# 3. Read lists. --tag-supplementary gives supplementary alignments the haplotype of their read.
for sample in tumor normal; do
  list=$OUT/haplotag/$sample.haplotags.tsv
  [[ -s $list ]] && continue
  timed whatshap_haplotag.$sample $WHATSHAP haplotag --reference $REF --ignore-read-groups --tag-supplementary \
    --output /dev/null --output-haplotag-list $list.part $PHASED $OUT/align/$sample.raw.bam \
    > $OUT/logs/whatshap_haplotag.$sample.log 2>&1
  mv $list.part $list
done

# 4. HP and PS tags into the raw and the error-corrected alignments.
for reads in raw ec; do
  for sample in tumor normal; do
    in=$OUT/align/$sample.$reads.bam
    hp=$OUT/align/$sample.$reads.hp.bam
    [[ -s $hp.bai ]] && continue
    timed tag.$sample.$reads $PYTHON $TAG_BAM $in $OUT/haplotag/$sample.haplotags.tsv $hp.part --threads 16 --with-ps \
      2> $OUT/logs/tag.$sample.$reads.log
    mv $hp.part $hp
    timed index.$sample.$reads.hp $SAMTOOLS index -@ 16 $hp
    rm $in $in.bai
  done
done
