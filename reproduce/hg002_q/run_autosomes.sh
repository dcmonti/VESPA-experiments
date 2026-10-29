#!/usr/bin/env bash
# HG002, autosomes: the chr1 analysis (run_all.sh) extended to chr1-chr22 (prepared and launched 2026-10-16).
# chrX and chrY are left out: HG002 is male, and the two-homolog model does not apply to them.
# The chr1 results of run_all.sh are not touched: alignments and summaries carry the tag "autosomes", and the
# per-chromosome tables go to $OUT/autosomes/.
# Usage: bash run_autosomes.sh
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${OUT:-$HG002/q_visibility}
CHROMS=$(echo chr{1..22})
export OUT CHROMS TAG=autosomes THREADS=48 JOBS=11
bash $HERE/01_data_align.sh
bash $HERE/02_haplotag.sh
printf "%s\n" $CHROMS | xargs -P 11 -I{} bash -c "$PYTHON $HERE/03_visibility.py $OUT {} autosomes 2> $OUT/logs/visibility.autosomes.{}.log"
printf "%s\n" $CHROMS | xargs -P 11 -I{} bash -c "$PYTHON $HERE/06_visibility_reads.py $OUT {} autosomes > $OUT/logs/visibility_reads.autosomes.{}.log 2>&1"
$PYTHON $HERE/07_summary.py $OUT/autosomes $CHROMS
bash $HERE/05_alpha.sh
