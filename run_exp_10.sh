#!/usr/bin/env bash
# Experiment 10. HG002, per-read detection rate q of LOWPOWER, autosomes.
# Germline heterozygous SVs of the GIAB v5.0q truth, PacBio Revio reads tagged by homolog with the truth phasing: per
# SV, the spanning reads of each homolog and those whose alignment shows the event; q and its checks.
# Inputs: reproduce/hg002_q/01_data_align.sh. Output: $HG002/q_visibility/autosomes.
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/config.sh"; source $LIB/common.sh

need $HG008/reads_ec_align/ref/GRCh38_GIABv3.chr1-22X.fa "run_exp_1.sh"   # the reads are aligned to the reference of the HG008 run
log "experiment 10: HG002 q"
bash $REPRO/hg002_q/run_autosomes.sh
log "experiment 10: done"
