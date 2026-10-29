#!/usr/bin/env bash
# Environment with PBSIM3 and ccs (once only).
set -euo pipefail
source "$(dirname "$0")/config.sh"
[[ -x $SIMENV/bin/pbsim && -x $SIMENV/bin/ccs ]] && { log "environment already present"; exit 0; }
$CONDA create -y -p $SIMENV -c conda-forge -c bioconda pbsim3 pbccs samtools
{ $SIMENV/bin/pbsim 2>&1 | head -3; $SIMENV/bin/ccs --version 2>&1 | head -1; } || true
find $SIMENV -name "QSHMM-RSII.model" | head -1
