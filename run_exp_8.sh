#!/usr/bin/env bash
# Experiment 8. Error correction on simulated tumor/normal pairs.
# A 5 Mb region of chr7 from the two haplotypes of the HG008-N verkko 2.2 assembly, known somatic variants in two
# clones, HiFi reads simulated with PBSIM3 + ccs at two depths and three purities; the tumor reads error-corrected with
# the NB model at three false-negative rates, with plain hifiasm, or not at all; then index, search and calls.
# Needs: the HG008-N verkko 2.2 assembly and hg38 (reproduce/sim_ec/config.sh). Output: $HG008/sim_ec_v2_2026-09-29.
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/config.sh"; source $LIB/common.sh

log "experiment 8: EC simulation"
bash $REPRO/sim_ec/run_all.sh
bash $REPRO/sim_ec/run_raw_err.sh
log "experiment 8: done"
