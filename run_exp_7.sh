#!/usr/bin/env bash
# Experiment 7. VAF of the ground-truth events, from svviz2 on the tumor raw reads.
#   HG008     V0.5 events of the clonal and subclonal regions, reads of experiment 1
#   COLO829   truth set, reads aligned with -s50 (experiment 3)
# The events are split in N chunks run in parallel; the events left without a report are run again in a final pass.
# Needs: experiments 1 and 3, the conda environment svviz2 (envs/). Output: <dataset run>/svviz2_vaf/vaf_events.tsv.
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/config.sh"; source $LIB/common.sh
V=$REPRO/vaf
N=8   # chunks, as in make_input.py, make_missing.py and make_chunks.py

need $HG008/reads_ec_align/bench/clonal_and_subclonal/out/HG008_support.tsv "run_exp_1.sh"
need $COLO829/reads_ec_align_s50/bench/min50/out/COLO829_support.tsv "run_exp_3.sh"

log "experiment 7: HG008"
$PYTHON $V/hg008/make_input.py
seq 0 $((N - 1)) | xargs -P $N -I{} bash $V/hg008/run.sh {}
$PYTHON $V/hg008/make_missing.py
seq 0 $((N - 1)) | xargs -P $N -I{} bash $V/hg008/run_missing.sh {}

log "experiment 7: COLO829"
O=$COLO829/reads_ec_align_s50/svviz2_vaf
mkdir -p $O/reports $O/logs $O/chunks
bash $V/colo829/make_input.sh
bash $V/colo829/run.sh > $O/logs/run.log 2>&1 || log "COLO829 single run incomplete, the chunks below finish it"
$PYTHON $V/colo829/make_chunks.py
seq 0 $((N - 1)) | xargs -P $N -I{} bash $V/colo829/run_chunk.sh {}

log "experiment 7: VAF table"
$PYTHON $V/svviz2_vaf_table.py
log "experiment 7: done"
