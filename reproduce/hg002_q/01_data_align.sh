#!/usr/bin/env bash
# HG002, chr1: data for the per-read detection rate q of LOWPOWER (Supplementary Methods, germline filtering).
# q = probability that the alignment of a normal read from the homolog carrying a germline heterozygous SV shows it.
# HG002 is not one of the evaluation datasets: q measured here is a property of the reads and of the aligner.
#
#   1. truth: GIAB HG002 v5.0q (T2T-HG002-Q100 v1.1), the dipcall VCF on GRCh38 (SNVs and SVs, phased in the frame of
#      the two assembly haplotypes) and the benchmark regions of the SV and small-variant sets
#   2. reads: PacBio Revio 2022Q4 HG002-rep1 (one SMRT cell); only the primary records on $CHROMS of the GRCh38 alignment
#      distributed by PacBio are fetched (remote BAM + index), and only to get the reads of those chromosomes
#   3. the raw reads are aligned again with minimap2 -ax map-hifi at its defaults to the reference of our runs, as in
#      reproduce/hg008/01_reads_ec_align.sh (no error correction: q is a property of the HiFi reads)
#
# Output, in $OUT: truth/, reads/hg002.$TAG.fq.gz (removed after the alignment when TAG != chr1),
#   align/hg002.$TAG.raw.bam(.bai), logs/
# Usage: bash 01_data_align.sh            (chr1)
#        CHROMS="chr1 ... chr22" TAG=autosomes bash 01_data_align.sh
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"

OUT=${OUT:-$HG002/q_visibility}
REF=$HG008/reads_ec_align/ref/GRCh38_GIABv3.chr1-22X.fa # the reference of the HG008 run
GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q
PB=https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1
PB_BAM=$PB/analysis/HG002.m84011_220902_175841_s1.GRCh38.bam
CHROMS=${CHROMS:-chr1}
TAG=${TAG:-chr1}
THREADS=${THREADS:-24}

mkdir -p $OUT/truth $OUT/reads $OUT/align $OUT/logs

fetch() { # url -> $OUT/truth/<basename>
  local dst=$OUT/truth/$(basename $1)
  [[ -s $dst ]] && return 0
  curl -fsSL --retry 3 -o $dst.part $1
  mv $dst.part $dst
}

# 1. truth
fetch $GIAB/dipcall_output/GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.vcf.gz
fetch $GIAB/dipcall_output/GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.vcf.gz.tbi
fetch $GIAB/HG002_GRCh38_v5.0q_stvar.benchmark.bed
fetch $GIAB/HG002_GRCh38_v5.0q_smvar.benchmark.bed
fetch $GIAB/NIST_HG002_v5.0q_variant-benchmarksets_README.md
fetch $GIAB/checksum.md5

# 2. reads of $CHROMS (primary records only, so that each read is written once)
fq=$OUT/reads/hg002.$TAG.fq.gz
bam=$OUT/align/hg002.$TAG.raw.bam
# One chromosome at a time, each retried: a long remote stream can break (libcurl error 92, HTTP/2 framing layer,
# on 2026-10-16 after 37 GB), and a failed chromosome is fetched again from its start. The reads of chr1 from the
# chr1-only run are the same records and are reused.
fetch_chrom() { # chrom -> $OUT/reads/hg002.$TAG.<chrom>.fq.gz
  local c=$1 dst=$OUT/reads/hg002.$TAG.$1.fq.gz try
  [[ -s $dst ]] && return 0
  if [[ $c == chr1 && -s $OUT/reads/hg002.chr1.fq.gz && $TAG != chr1 ]]; then
    ln $OUT/reads/hg002.chr1.fq.gz $dst; return 0
  fi
  $SV_BIN/samtools view -H $PB_BAM | grep -q "SN:$c\b" || { echo "no $c in $PB_BAM" >&2; return 1; }
  for try in 1 2 3 4; do
    if $SV_BIN/samtools view -u -F 0x900 $PB_BAM $c 2>> $OUT/logs/fetch.$TAG.log \
         | $SV_BIN/samtools fastq -@ 4 - 2>> $OUT/logs/fetch.$TAG.log | gzip -1 > $dst.part; then
      mv $dst.part $dst; return 0
    fi
    echo "fetch of $c failed (try $try)" >> $OUT/logs/fetch.$TAG.log; rm -f $dst.part; sleep 60
  done
  return 1
}
if [[ ! -s $fq && ! -s $bam.bai ]]; then
  if [[ $TAG == chr1 ]]; then
    $SV_BIN/samtools view -H $PB_BAM | grep -q "SN:chr1\b" || { echo "no chr1 in $PB_BAM" >&2; exit 1; }
    $SV_BIN/samtools view -u -F 0x900 $PB_BAM chr1 2> $OUT/logs/fetch.log \
      | $SV_BIN/samtools fastq -@ 4 - 2>> $OUT/logs/fetch.log | gzip -1 > $fq.part
    mv $fq.part $fq
  else
    for c in $CHROMS; do fetch_chrom $c; done
    cat $(for c in $CHROMS; do echo $OUT/reads/hg002.$TAG.$c.fq.gz; done) > $fq.part
    mv $fq.part $fq
    rm -f $OUT/reads/hg002.$TAG.chr*.fq.gz
  fi
fi

# 3. alignment, as in our runs
if [[ ! -s $bam.bai ]]; then
  log=$OUT/logs/minimap2.log; [[ $TAG == chr1 ]] || log=$OUT/logs/minimap2.$TAG.log
  $SV_BIN/minimap2 -ax map-hifi -t $THREADS $REF $fq 2> $log \
    | $SV_BIN/samtools sort -@ 8 -m 4G -T $OUT/align/sort.$TAG -o $bam.part -
  mv $bam.part $bam
  $SV_BIN/samtools index -@ 8 $bam
  [[ $TAG == chr1 ]] || rm -f $fq      # the reads can be fetched again; chr1 keeps its file as before
fi
echo done
