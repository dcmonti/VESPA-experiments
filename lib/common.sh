# Helpers of the run_exp_N.sh scripts.

log() { echo "[$(date '+%F %T')] $*"; }

# need FILE EXPERIMENT: stop unless FILE, an output of EXPERIMENT, exists.
need() { [[ -s $1 ]] || { echo "missing $1: run $2 first" >&2; exit 1; }; }
