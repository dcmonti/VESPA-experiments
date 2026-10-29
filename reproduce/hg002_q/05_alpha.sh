#!/usr/bin/env bash
# HG002: over-dispersion alpha of the read depth, as for the two normals (reproduce/depth_dispersion), on the raw
# alignment of 01_data_align.sh. The positions are those of depth_dispersion (same seed, genome-wide) that fall on
# $CHROMS, so they are the same positions as in the HG008 and COLO829 estimates. alpha is a parameter of the error
# correction (negative binomial of the threshold T), not of LOWPOWER; HG002 gives a third, independent run.
#
# Output, in $OUT/depth_dispersion: positions.$TAG.tsv.gz, depth.$TAG.tsv.gz, summary.$TAG.tsv (chr1: summary.tsv)
# Usage: bash 05_alpha.sh      or      CHROMS="chr1 ... chr22" TAG=autosomes bash 05_alpha.sh
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"

OUT=${OUT:-$HG002/q_visibility}
DD=$(cd "$(dirname "$0")/../depth_dispersion" && pwd)
REF=$REFS/GRCh38_GIABv3.fasta.gz
EXCL=$REFS/grch38_centromers.bed
TRF=$REFS/human_GRCh38_no_alt_analysis_set.trf.bed
CHROMS=${CHROMS:-chr1}
TAG=${TAG:-chr1}

D=$OUT/depth_dispersion
mkdir -p $D
if [[ ! -s $D/positions.$TAG.tsv.gz ]]; then
  [[ -s $D/positions.all.tsv.gz ]] || \
    $PYTHON $DD/dispersion.py positions --ref $REF --exclude $EXCL --n 100000 --out $D/positions.all.tsv.gz
  zcat $D/positions.all.tsv.gz | awk -v c="$CHROMS" 'BEGIN{n = split(c, a, " "); for (i = 1; i <= n; i++) k[a[i]] = 1}
                                    NR == 1 || ($1 in k)' | gzip > $D/positions.$TAG.tsv.gz
fi
[[ -s $D/depth.$TAG.tsv.gz ]] || \
  $PYTHON $DD/dispersion.py depth --bam $OUT/align/hg002.$TAG.raw.bam --positions $D/positions.$TAG.tsv.gz \
    --threads 16 --out $D/depth.$TAG.tsv.gz
sum=$D/summary.$TAG.tsv; [[ $TAG == chr1 ]] && sum=$D/summary.tsv
$PYTHON $DD/dispersion.py fit --depth $D/depth.$TAG.tsv.gz --trf $TRF --sample hg002.$TAG --boot 200 --out $sum
cat $sum
