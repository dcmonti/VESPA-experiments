#!/usr/bin/env bash
# Experiment 4. Effect of each germline filter.
# VESPA without the normal, then each germline filter added in turn, benchmarked with Minda: HG008 at s = 5,
# COLO829 at s = 8 with the reads aligned with -s50.
# Needs: experiments 1 and 3. Output: <dataset run>/germline_filters/summary.tsv.
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/config.sh"; source $LIB/common.sh

need $HG008/reads_ec_align/bench/clonal_and_subclonal/summary.tsv "run_exp_1.sh"
need $COLO829/reads_ec_align_s50/bench/min50/summary.tsv "run_exp_3.sh"
log "experiment 4: germline filters"
bash $REPRO/germline_filters/run.sh hg008 colo829
log "experiment 4: done"
