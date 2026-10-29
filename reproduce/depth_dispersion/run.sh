#!/usr/bin/env bash
# Over-dispersion alpha of the read depth of the two normals (raw alignments of reads_ec_align), for the negative
# binomial of the error-correction threshold. Same random autosomal positions for both samples (fixed seed).
# Method in dispersion.py.
#
# Output, in <dataset>/reads_ec_align/depth_dispersion/:
#   positions.tsv.gz      sampled positions (0-based)
#   normal.depth.tsv.gz   depth of the raw normal at each position
#   summary.tsv           mean, variance, alpha (moments and ML) with bootstrap 95% intervals,
#                         on all positions and outside tandem repeats
#
# Usage: bash run.sh
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"
HERE=$(cd "$(dirname "$0")" && pwd)
REF=$REFS/GRCh38_GIABv3.fasta.gz
EXCL=$REFS/grch38_centromers.bed
TRF=$REFS/human_GRCh38_no_alt_analysis_set.trf.bed

declare -A RUN=(
  [hg008]=$HG008/reads_ec_align
  [colo829]=$COLO829/reads_ec_align
)

for ds in hg008 colo829; do
  OUT=${RUN[$ds]}/depth_dispersion
  mkdir -p $OUT
  [[ -s $OUT/positions.tsv.gz ]] || \
    $PYTHON $HERE/dispersion.py positions --ref $REF --exclude $EXCL --n 100000 --out $OUT/positions.tsv.gz
  [[ -s $OUT/normal.depth.tsv.gz ]] || \
    /usr/bin/time -v -o $OUT/time.depth.txt \
    $PYTHON $HERE/dispersion.py depth --bam ${RUN[$ds]}/align/normal.raw.hp.bam \
      --positions $OUT/positions.tsv.gz --threads 22 --out $OUT/normal.depth.tsv.gz
  $PYTHON $HERE/dispersion.py fit --depth $OUT/normal.depth.tsv.gz --trf $TRF --sample $ds.normal \
    --boot 200 --out $OUT/summary.tsv
  cat $OUT/summary.tsv
done
