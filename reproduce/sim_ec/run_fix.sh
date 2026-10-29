#!/usr/bin/env bash
# Fix pass over all cells (08_retag_recall.sh), then the merged tables in table/ (small, bp, calls,
# noise, attr, attr_background).
# USAGE: PAR=3 nohup bash run_fix.sh > logs/run_fix.log 2>&1 &
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"; source ./config.sh; L=$D/logs
log "start of fix: $EC_SETS, $PAR in parallel"
for s in $(samples | grep -v normal); do for e in $EC_SETS; do echo "$s $e"; done; done \
  | xargs -P $PAR -L1 bash -c "bash 08_retag_recall.sh \$0 \$1 > $L/08_fix.\$0.\$1.log 2>&1 || echo \"cell \$0 \$1 FAILED\""
for t in small bp calls noise attr attr_background; do
  f=$(ls $D/cells/*/*/eval/$t.tsv | head -1); head -1 $f > $D/table/$t.tsv
  for g in $D/cells/*/*/eval/$t.tsv; do tail -n +2 $g; done >> $D/table/$t.tsv
done
log "done: tables in $D/table"
