#!/usr/bin/env bash
# Every experiment, in an order that satisfies the dependencies (see each run_exp_N.sh). Each step skips the
# outputs it already finds, so after a failure the same command resumes where it stopped.
# Usage: nohup bash run_exp_all.sh > run_exp_all.log 2>&1 &
set -euo pipefail
HERE=$(dirname "$(readlink -f "$0")")
source $HERE/config.sh; source $LIB/common.sh

for n in 1 2 3 4 5 6 7 8 9 10; do
  bash $HERE/run_exp_$n.sh
done
log "all experiments done"
