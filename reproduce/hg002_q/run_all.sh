#!/usr/bin/env bash
# HG002, chr1: per-read detection rate q of LOWPOWER and over-dispersion alpha of the depth, on a dataset that is not
# one of the evaluation datasets (prepared 2026-10-16; see the header of each step).
#   01_data_align.sh  truth (GIAB v5.0q, T2T-HG002-Q100 v1.1) and the chr1 reads of PacBio Revio HG002-rep1, minimap2
#   02_haplotag.sh    reads tagged by homolog with the truth phasing
#   03_visibility.py  per SV, spanning reads of each homolog and those that show it (strict / loose criterion)
#   04_fit.py         q, its checks (independence, other homolog, homozygous SVs)
#   05_alpha.sh       alpha of the depth on chr1
# Usage: bash run_all.sh
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${OUT:-$HG002/q_visibility}
export OUT
bash $HERE/01_data_align.sh
bash $HERE/02_haplotag.sh
$PYTHON $HERE/03_visibility.py $OUT
$PYTHON $HERE/04_fit.py $OUT
bash $HERE/05_alpha.sh
