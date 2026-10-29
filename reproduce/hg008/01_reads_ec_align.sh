#!/usr/bin/env bash
# HG008: read extraction, per-chromosome error correction, alignment.
#
# Input: the GIAB HG008 PacBio Revio BAMs aligned by NIST to GRCh38-GIABv3
# (pbmm2 1.10.0, --preset HIFI --unmapped), from
# https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/PacBio_Revio_20240125/
#   tumor  HG008-T_PacBio-HiFi-Revio_20240125_116x_GRCh38-GIABv3.bam    md5 60fe2f126ecfc1e1085343d744e57ad0
#   normal HG008-N-P_PacBio-HiFi-Revio_20240125_35x_GRCh38-GIABv3.bam   md5 83b03d37f23db9405246aedfedc1c0da
# Reference: GRCh38-GIABv3, the reference of the NIST alignments and of the HG008 benchmark, from
# https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/
#   GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz   md5 939ce19062d1462c09b88c55faca4d76
# Chromosomes: chr1-22 and chrX. chrY is excluded because the HG008 donor is female.
#
# Steps:
#   1. restrict the reference to the selected chromosomes
#   2. extract the reads whose primary alignment is on a selected chromosome, one FASTQ per chromosome
#   3. error-correct each chromosome separately with hifiasm, then merge the chromosomes
#   4. align the raw reads with minimap2 -ax map-hifi (default parameters)
#   5. align the error-corrected reads in the same way
#
# Output, in $OUT:
#   ref/GRCh38_GIABv3.chr1-22X.fa          reference used for the alignments
#   ec/{tumor,normal}.ec.fa.gz             error-corrected reads
#   ec/normal.r_utg.fa.gz                  unitigs of the normal (hifiasm r_utg graph), input of the VESPA index
#   align/{tumor,normal}.{raw,ec}.bam      alignments, sorted and indexed
#   time/<step>.txt                        GNU time (/usr/bin/time -v) report of each step
# The raw reads extracted in step 2 are deleted once aligned.
#
# Usage: bash 01_reads_ec_align.sh
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"

TUMOR_BAM=$HG008/jan2024_2026-09-20/raw/HG008-T.jan2024.bam
NORMAL_BAM=$HG008/jan2024_2026-09-20/raw/HG008-N-P.jan2024.bam
GIABV3=$REFS/GRCh38_GIABv3.fasta.gz
OUT=${OUT:-$HG008/reads_ec_align}

THREADS=48
CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"
TUMOR_HET_COV=6
NORMAL_HET_COV=4

mkdir -p $OUT/ref $OUT/reads $OUT/ec $OUT/align $OUT/logs $OUT/time
REF=$OUT/ref/GRCh38_GIABv3.chr1-22X.fa

# timed NAME COMMAND...: run COMMAND under GNU time, report in $OUT/time/NAME.txt. A pipeline is timed as one
# step by passing it to bash -o pipefail -c.
timed() { local name=$1; shift; /usr/bin/time -v -o $OUT/time/$name.txt "$@"; }

input_bam() { if [[ $1 == tumor ]]; then echo $TUMOR_BAM; else echo $NORMAL_BAM; fi; }
het_cov()   { if [[ $1 == tumor ]]; then echo $TUMOR_HET_COV; else echo $NORMAL_HET_COV; fi; }

# 1. Reference restricted to the selected chromosomes, in the order of $CHROMS.
if [[ ! -s $REF.fai ]]; then
  timed ref_extract $SAMTOOLS faidx $GIABV3 $CHROMS > $REF.part
  mv $REF.part $REF
  timed ref_index $SAMTOOLS faidx $REF
fi

# 2. Extract reads: -F 0x900 drops secondary and supplementary records, so each read is taken once,
#    on the chromosome of its primary alignment. Unmapped reads and reads on other contigs are not extracted.
for sample in tumor normal; do
  [[ -s $OUT/align/$sample.raw.bam.bai ]] && continue   # already aligned, reads deleted
  for chrom in $CHROMS; do
    fq=$OUT/reads/$sample.$chrom.fq.gz
    [[ -s $fq ]] && continue
    timed extract.$sample.$chrom bash -o pipefail -c "$SAMTOOLS view -u -F 0x900 $(input_bam $sample) $chrom | $SAMTOOLS fastq - | gzip -1 > $fq.part"
    mv $fq.part $fq
  done
done

# 3. Error correction, one chromosome at a time, then one file per sample.
#    Tumor: --bin-only, only the corrected reads are needed.
#    Normal: full run, which also writes the unitig graph; its unitigs are kept as FASTA, each name
#    prefixed with the chromosome because hifiasm gives the same names in every run.
for sample in tumor normal; do
  ec=$OUT/ec/$sample.ec.fa.gz
  [[ -s $ec ]] && continue
  for chrom in $CHROMS; do
    fq=$OUT/reads/$sample.$chrom.fq.gz
    ec_chrom=$OUT/ec/$sample.$chrom.ec.fa.gz
    work=$OUT/ec/$sample.$chrom.work
    [[ -s $ec_chrom ]] && continue
    bin_only=""; [[ $sample == tumor ]] && bin_only=--bin-only
    rm -rf $work; mkdir -p $work
    timed hifiasm.$sample.$chrom $HIFIASM -t $THREADS --write-ec --het-cov $(het_cov $sample) $bin_only -o $work/ec $fq \
      > $OUT/logs/hifiasm.$sample.$chrom.log 2>&1
    if [[ $sample == normal ]]; then
      awk -v p=${chrom}_ '/^S/ {print ">" p $2; print $3}' $work/ec.bp.r_utg.gfa | gzip -1 > $OUT/ec/normal.$chrom.r_utg.fa.gz
    fi
    gzip -1 -c $work/ec.ec.fa > $ec_chrom.part
    mv $ec_chrom.part $ec_chrom
    rm -rf $work
  done
  # concatenated gzip files are a valid gzip file
  if [[ $sample == normal ]]; then
    cat $(for c in $CHROMS; do echo $OUT/ec/normal.$c.r_utg.fa.gz; done) > $OUT/ec/normal.r_utg.fa.gz
  fi
  cat $(for c in $CHROMS; do echo $OUT/ec/$sample.$c.ec.fa.gz; done) > $ec.part
  mv $ec.part $ec
  rm -f $OUT/ec/$sample.chr*.ec.fa.gz $OUT/ec/$sample.chr*.r_utg.fa.gz
done

# 4-5. Alignment with minimap2 default HiFi preset, raw reads then error-corrected reads.
for reads in raw ec; do
  for sample in tumor normal; do
    bam=$OUT/align/$sample.$reads.bam
    [[ -s $bam.bai ]] && continue
    # one line, space-separated: the list goes inside the bash -c command string
    if [[ $reads == raw ]]; then files=$(printf "$OUT/reads/$sample.%s.fq.gz " $CHROMS)
    else                         files=$OUT/ec/$sample.ec.fa.gz; fi
    timed align.$sample.$reads bash -o pipefail -c "zcat $files \
      | $MINIMAP2 -ax map-hifi -t $THREADS $REF - 2> $OUT/logs/minimap2.$sample.$reads.log \
      | $SAMTOOLS sort -@ 16 -m 4G -T $OUT/align/$sample.$reads.sort -o $bam.part -"
    mv $bam.part $bam
    timed index.$sample.$reads $SAMTOOLS index -@ 16 $bam
    if [[ $reads == raw ]]; then rm $files; fi
  done
done
