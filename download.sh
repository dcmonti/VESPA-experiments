#!/usr/bin/env bash
# Inputs of the experiments, downloaded to the paths the scripts read (config.sh), each checked against its md5.
# Published call sets are not downloaded: see README.md. HG002 (experiment 10) downloads its own inputs.
#   refs     GRCh38-GIABv3, tandem repeats of Severus, centromeres and telomeres                  (all experiments)
#   hg008    HG008 PacBio Revio reads, V0.5 SV and V0.3 small-variant benchmarks, HG008-N v6.3        (1, 4-7)
#   colo829  COLO829 PacBio Revio reads, truth set                                                      (2-4, 7, 9)
#   sim      hg38 primary chromosomes and the HG008-N verkko 2.2 haplotypes                             (8)
# Usage: bash download.sh [refs] [hg008] [colo829] [sim]      (no argument: all)
#        CHECK=1 bash download.sh ...                       (only check that every URL answers, download nothing)
set -euo pipefail
HERE=$(dirname "$(readlink -f "$0")")
source $HERE/config.sh; source $LIB/common.sh

GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab
HG008_FTP=$GIAB/data_somatic/HG008/Liss_lab
NIST_S3=https://nist-giab.s3.us-east-1.amazonaws.com/giab_tumor-normal/analysis/HG008/NIST_asm_dev
PACBIO=https://downloads.pacbcloud.com/public/revio/2023Q2/COLO829
ZENODO=https://zenodo.org/records
CHECK=${CHECK:-0}

# fetch URL DEST [MD5]: download unless DEST exists, then check the md5 when one is given.
fetch() {
  local url=$1 dst=$2 md5=${3:-}
  if [[ $CHECK == 1 ]]; then
    echo -e "$(curl -sIL -o /dev/null -w '%{http_code}' "$url")\t$url"; return 0
  fi
  mkdir -p "$(dirname $dst)"
  if [[ ! -s $dst ]]; then
    log "download $url"
    curl -fL --retry 5 -C - -o $dst.part "$url"
    mv $dst.part $dst
  fi
  [[ -z $md5 ]] || echo "$md5  $dst" | md5sum -c --quiet - || { echo "md5 mismatch: $dst" >&2; exit 1; }
}

refs() {
  local f=GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
  fetch $GIAB/release/references/GRCh38/$f     $REFS/GRCh38_GIABv3.fasta.gz 939ce19062d1462c09b88c55faca4d76
  fetch $GIAB/release/references/GRCh38/$f.fai $REFS/GRCh38_GIABv3.fasta.gz.fai
  fetch $GIAB/release/references/GRCh38/$f.gzi $REFS/GRCh38_GIABv3.fasta.gz.gzi
  # Severus 1.7, vntrs/: the file used by Minda and given to Severus and Sniffles2
  fetch https://raw.githubusercontent.com/KolmogorovLab/Severus/1.7/vntrs/human_GRCh38_no_alt_analysis_set.trf.bed \
        $REFS/human_GRCh38_no_alt_analysis_set.trf.bed 86d75a85add01f940c4d5abc4bd596b9
  # no public source: distributed with this repository
  [[ $CHECK == 1 || -s $REFS/grch38_centromers.bed ]] || cp $HERE/data/grch38_centromers.bed $REFS/
}

hg008() {
  local raw=$HG008/jan2024_2026-09-20/raw r=$HG008_FTP/PacBio_Revio_20240125
  fetch $r/HG008-T_PacBio-HiFi-Revio_20240125_116x_GRCh38-GIABv3.bam       $raw/HG008-T.jan2024.bam   60fe2f126ecfc1e1085343d744e57ad0
  fetch $r/HG008-T_PacBio-HiFi-Revio_20240125_116x_GRCh38-GIABv3.bam.bai   $raw/HG008-T.jan2024.bam.bai
  fetch $r/HG008-N-P_PacBio-HiFi-Revio_20240125_35x_GRCh38-GIABv3.bam      $raw/HG008-N-P.jan2024.bam 83b03d37f23db9405246aedfedc1c0da
  fetch $r/HG008-N-P_PacBio-HiFi-Revio_20240125_35x_GRCh38-GIABv3.bam.bai  $raw/HG008-N-P.jan2024.bam.bai

  # V0.5 somatic SV draft benchmark; the VCF is used uncompressed
  local v=$HG008_FTP/analysis/NIST_HG008-T_somatic-stvar-CNV_DraftBenchmark_V0.5-20260318 g=$HG008/ground_truths/v0.5
  local vcf=GRCh38_HG008-T-V0.5_somatic-stvar_ALL.draftbenchmark.vcf
  fetch $v/$vcf.gz $g/$vcf.gz
  if [[ $CHECK == 0 ]]; then
    [[ -s $g/$vcf ]] || zcat $g/$vcf.gz > $g/$vcf
    echo "778225b2dff9a74797c6632d03113544  $g/$vcf" | md5sum -c --quiet -
  fi
  fetch $v/GRCh38_HG008-T-V0.5_somatic-stvar-clonal.draftbenchmark.bed             $g/GRCh38_HG008-T-V0.5_somatic-stvar-clonal.draftbenchmark.bed
  fetch $v/GRCh38_HG008-T-V0.5_somatic-stvar-clonal_and_subclonal.draftbenchmark.bed $g/GRCh38_HG008-T-V0.5_somatic-stvar-clonal_and_subclonal.draftbenchmark.bed
  fetch $v/README.md $g/README.md

  # V0.3 somatic small-variant draft benchmark (experiment 6)
  local s=$HG008_FTP/analysis/NIST_HG008-T_somatic-smvar_DraftBenchmark_V0.3-20260425 m=$HG008/ground_truths/smvar_v0.3
  for f in HG008-T_somatic_smvar_benchmark_v0.3_all.bed HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz \
           HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz.tbi README.md; do
    fetch $s/$f $m/$f
  done

  # HG008-N curated assembly v6.3 (experiment 6)
  fetch $NIST_S3/HG008N_v6.3/HG008N_v6.3.fasta.gz $HG008/refs/HG008N_v6.3/HG008N_v6.3.fasta.gz a5de2e51f0ad423d236036ed61e2101c
  fetch $NIST_S3/HG008N_v6.3/HG008N_v6.3.gaps.bed $HG008/refs/HG008N_v6.3/HG008N_v6.3.gaps.bed 547379a0483e3d807b13f11479fc989a
}

# HG008-N verkko 2.2, haplotypes 1 and 2 (experiment 8); diploid.fa.gz: both, names prefixed hap1_ / hap2_
verkko22() {
  local k=$HG008_FTP/analysis/Verkko_assemblies_05162024/HG008_N_asm_hifiherrohic_verkko2.2_20250218 d=$HG008/refs/HG008N_verkko2.2_20250218
  fetch $k/HG008N_verkko-assembly_20250218.haplotype1.fasta.gz $d/HG008N_verkko-assembly_20250218.haplotype1.fasta.gz 446b7d66fb7042c44e21d095c69fd8be
  fetch $k/HG008N_verkko-assembly_20250218.haplotype2.fasta.gz $d/HG008N_verkko-assembly_20250218.haplotype2.fasta.gz cf8524a5836168f0758cb8155889d034
  if [[ $CHECK == 0 && ! -s $d/diploid.fa.gz.fai ]]; then
    { zcat $d/HG008N_verkko-assembly_20250218.haplotype1.fasta.gz | sed 's/^>/>hap1_/'
      zcat $d/HG008N_verkko-assembly_20250218.haplotype2.fasta.gz | sed 's/^>/>hap2_/'; } | $SV_BIN/bgzip -@ 8 > $d/diploid.fa.gz
    $SAMTOOLS faidx $d/diploid.fa.gz
  fi
}

colo829() {
  fetch $PACBIO/COLO829/analysis/COLO829.GRCh38.bam           $COLO829/samples/tumor.colo829.bam  ccbcffc443f6417d4aabe6c5eb944211
  fetch $PACBIO/COLO829/analysis/COLO829.GRCh38.bam.bai       $COLO829/samples/tumor.colo829.bam.bai
  fetch $PACBIO/COLO829-BL/analysis/COLO829-BL.GRCh38.bam     $COLO829/samples/normal.colo829.bam c7c234144688b45d781940e13069d817
  fetch $PACBIO/COLO829-BL/analysis/COLO829-BL.GRCh38.bam.bai $COLO829/samples/normal.colo829.bam.bai
  # truth set of Valle-Inclan et al. 2022 on GRCh38; breakend pairs collapsed into symbolic records for truvari
  fetch $ZENODO/7515830/files/truthset_somaticSVs_COLO829_hg38lifted.vcf $COLO829/gt/truthset_COLO829_hg38lifted.vcf 504e4ca4bd285706c261f93be1f33a95
  [[ $CHECK == 1 || -s $COLO829/gt/truthset_COLO829_hg38lifted.symbolic.vcf ]] || \
    $PYTHON $LIB/bnd_pairs_to_symbolic.py $COLO829/gt/truthset_COLO829_hg38lifted.vcf -o $COLO829/gt/truthset_COLO829_hg38lifted.symbolic.vcf
}

sim() {
  # UCSC hg38, chr1-22, X, Y, M
  local h=$HG008/refs
  fetch https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz $h/hg38.fa.gz
  if [[ $CHECK == 0 && ! -s $h/hg38.primary.fa.gz.fai ]]; then
    zcat $h/hg38.fa.gz | $SV_BIN/bgzip -@ 8 > $h/hg38.fa.bgz && $SAMTOOLS faidx $h/hg38.fa.bgz
    $SAMTOOLS faidx $h/hg38.fa.bgz $(echo chr{1..22} chrX chrY chrM | tr ' ' '\n' | sort) | $SV_BIN/bgzip -@ 8 > $h/hg38.primary.fa.gz
    $SAMTOOLS faidx $h/hg38.primary.fa.gz
    rm -f $h/hg38.fa.bgz $h/hg38.fa.bgz.fai $h/hg38.fa.bgz.gzi
  fi
  verkko22
}

for target in ${@:-refs hg008 colo829 sim}; do
  log "inputs: $target"; $target
done
log "inputs: done"
