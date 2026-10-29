#!/usr/bin/env bash
# COLO829, second alignment (minimap2 -s50, see s50_align_tag.sh): alignment and tag transfer, then the calls and
# the benchmark of the first run, unchanged, with OUT pointing to the new directory. Each step skips the outputs it
# already finds, so after a failure the same command resumes from the step that failed.
# Run it after run_all.sh has finished (it reads the BAMs, read lists, SFS and reference of the first run).
#
# Usage: nohup bash s50_run_all.sh > $COLO829/reads_ec_align_s50.run_all.log 2>&1 &
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"

export OUT_DEFAULT=${OUT_DEFAULT:-$COLO829/reads_ec_align}
export OUT=${OUT:-$COLO829/reads_ec_align_s50}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p $OUT/logs

for step in s50_align_tag 04_call 05_bench; do
  echo "[$(date '+%F %T')] $step"
  bash $HERE/$step.sh > $OUT/logs/$step.log 2>&1
done
echo "[$(date '+%F %T')] done"
