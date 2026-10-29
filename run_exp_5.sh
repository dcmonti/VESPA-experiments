#!/usr/bin/env bash
# Experiment 5. HG008 with and without homopolymer compression (HPC) and the 500 bp cap on SFS length.
# Everything else as experiment 1, calls at s = 5, benchmark with Minda.
# Needs: experiment 1. Output: $HG008/reads_ec_align/{nohpc,sfs_lencap}.
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/config.sh"; source $LIB/common.sh

need $HG008/reads_ec_align/bench/clonal_and_subclonal/summary.tsv "run_exp_1.sh"
log "experiment 5: no HPC";               bash $REPRO/hg008_nohpc/run.sh
log "experiment 5: no HPC, no length cap"; bash $REPRO/hg008_nohpc/run_nocap.sh
log "experiment 5: no length cap";        bash $REPRO/hg008_sfs_lencap/run.sh
log "experiment 5: done"
