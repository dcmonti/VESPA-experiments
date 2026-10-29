#!/usr/bin/env bash
# HG008: all the steps, one after the other. Each step skips the outputs it already finds in $OUT, so after a
# failure the same command resumes from the step that failed.
#
# Usage: nohup bash run_all.sh > run_all.log 2>&1 &
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"

export OUT=${OUT:-$HG008/reads_ec_align}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p $OUT/logs

for step in 01_reads_ec_align 02_haplotag 03_vespa_index_search 04_call 05_bench; do
  echo "[$(date '+%F %T')] $step"
  bash $HERE/$step.sh > $OUT/logs/$step.log 2>&1
done
echo "[$(date '+%F %T')] done"
