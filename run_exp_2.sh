#!/usr/bin/env bash
# Experiment 2. COLO829, comparison of the callers, reads aligned with the minimap2 default (-s200).
# The steps of experiment 1 on COLO829, benchmarked genome-wide against the truth set of Valle-Inclan et al.,
# together with the published callsets of Qin et al. and Keskus et al.
# Inputs: reproduce/colo829/01_reads_ec_align.sh. Output: $COLO829/reads_ec_align, bench/min50/summary.tsv.
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/config.sh"; source $LIB/common.sh

log "experiment 2: COLO829 -s200"
bash $REPRO/colo829/run_all.sh
log "experiment 2: done"
