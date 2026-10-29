#!/usr/bin/env bash
# Addition of 2026-09-30 to the v2 simulation: (1) raw cells, tumor without EC, control for the whole
# row (SFS, cluster, attribution); (2) residual error per kb relative to the sequence of origin (08_err.py) on
# all cells, raw included; then the TSVs of table/ and the tables. The cells with EC are not rerun.
# USAGE: THREADS=16 PAR=3 nohup bash run_raw_err.sh > logs/run_raw_err.log 2>&1 &
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"; source ./config.sh; L=$D/logs
log "start: raw cells and residual error in $D, $PAR in parallel, $THREADS threads"
samples | grep -v normal | xargs -P $PAR -I{} bash -c "bash 05_cell.sh {} raw > $L/05_cell.{}.raw.log 2>&1 \
  && python3 06_eval.py $D {} raw '$WEIGHTS' > $L/06_eval.{}.raw.log 2>&1 && python3 07_attrib.py $D {} raw > $L/07_attrib.{}.raw.log 2>&1 \
  || echo 'cell {} raw FAILED'"
log "raw cells done, residual error"
for s in $(samples | grep -v normal); do for e in $EC_SETS; do echo "$s $e"; done; done \
  | xargs -P $PAR -L1 bash -c "[[ -s $D/cells/\$0/\$1/eval/err.tsv ]] || python3 08_err.py $D \$0 \$1 $THREADS > $L/08_err.\$0.\$1.log 2>&1 || echo \"error \$0 \$1 FAILED\""
mkdir -p $D/table
for t in small bp calls noise attr attr_background err; do
  f=$(ls $D/cells/*/*/eval/$t.tsv | head -1); head -1 $f > $D/table/$t.tsv
  for g in $D/cells/*/*/eval/$t.tsv; do tail -n +2 $g; done >> $D/table/$t.tsv
done
log "done: $D/table/"
