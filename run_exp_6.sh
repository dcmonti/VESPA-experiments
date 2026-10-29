#!/usr/bin/env bash
# Experiment 6. HG008, somatic signal discovery with SFS.
#   bp_coverage/        VESPA call without and with the normal, --call-provenance: the stage each ground-truth
#                       breakpoint reaches
#   germline_cascade/   indexes built from growing germline sets (reference only, + unitigs) and from the reference
#                       + the curated HG008-N v6.3 assembly; SFS, calls and Minda for each
#   sfs_signal/         per-breakpoint SFS coverage, classes of the SFS, the breakpoints no SFS covers
# Needs: experiment 1. Output: $HG008/reads_ec_align/{bp_coverage,germline_cascade,sfs_signal}.
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/config.sh"; source $LIB/common.sh

need $HG008/reads_ec_align/bench/clonal_and_subclonal/summary.tsv "run_exp_1.sh"
log "experiment 6: breakpoint coverage"; bash $REPRO/bp_coverage/run.sh
log "experiment 6: germline cascade";    bash $REPRO/germline_cascade/run_all.sh
log "experiment 6: HG008-N v6.3";      bash $REPRO/germline_cascade/06_v63.sh
# in order of dependence: bp_points -> lost_points -> tr_alleles -> tr_sfs_somatic; sfs_classes -> smvar_pileup ->
# smvar_join, cluster_to_sv, hp_runs
for s in bp_points lost_points tr_alleles tr_sfs_somatic sfs_classes smvar_pileup smvar_join cluster_to_sv hp_runs trf_fraction; do
  log "experiment 6: $s"; $PYTHON $REPRO/sfs_signal/$s.py
done
log "experiment 6: done"
