#!/usr/bin/env bash
# Experiment 1. HG008, comparison of the callers.
# Error correction, alignment, haplotagging, index and SFS search, calls of VESPA, Severus, SAVANA, nanomonsv and
# Sniffles2 (default run and threshold grid), benchmark with Minda and truvari against the V0.5 draft benchmark in
# the clonal and in the clonal and subclonal regions, together with the published HiFi callsets of the Minda ensemble.
# Inputs: reproduce/hg008/01_reads_ec_align.sh. Output: $HG008/reads_ec_align, bench/<regions>/summary.tsv.
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/config.sh"; source $LIB/common.sh

log "experiment 1: HG008"
bash $REPRO/hg008/run_all.sh
log "experiment 1: done"
