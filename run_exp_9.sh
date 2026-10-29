#!/usr/bin/env bash
# Experiment 9. Over-dispersion alpha of the read depth of the normal.
# Depth of the raw normal alignments at 100 000 random autosomal positions outside N bases, centromeres and
# telomeres, negative binomial fit, HG008 and COLO829.
# Needs: experiments 1 and 2. Output: <dataset run>/depth_dispersion/summary.tsv.
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/config.sh"; source $LIB/common.sh

need $HG008/reads_ec_align/align/normal.raw.hp.bam "run_exp_1.sh"
need $COLO829/reads_ec_align/align/normal.raw.hp.bam "run_exp_2.sh"
log "experiment 9: depth over-dispersion"
bash $REPRO/depth_dispersion/run.sh
log "experiment 9: done"
