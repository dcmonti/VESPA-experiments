#!/usr/bin/env bash
# Experiment 3. COLO829, comparison of the callers, reads aligned with minimap2 -s50.
# Raw and error-corrected reads of experiment 2 realigned with -s50, haplotags transferred, then the calls and the
# benchmark of experiment 2 unchanged.
# Needs: experiment 2. Output: $COLO829/reads_ec_align_s50, bench/min50/summary.tsv.
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/config.sh"; source $LIB/common.sh

need $COLO829/reads_ec_align/bench/min50/summary.tsv "run_exp_2.sh"
log "experiment 3: COLO829 -s50"
bash $REPRO/colo829/s50_run_all.sh
log "experiment 3: done"
