#!/usr/bin/env bash
# Full simulation chain. Each step skips what already exists. The outputs go to $D
# (config.sh; SIM_DIR for a different folder, for example a test run).
# USAGE: THREADS=16 PAR=3 nohup bash run_all.sh > logs/run_all.log 2>&1 &
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"; source ./config.sh; mkdir -p $D/logs; L=$D/logs
log "start: $D, region $REGION, depths $DEPTHS, purities $PURITIES, EC $EC_SETS, $PAR in parallel, $THREADS threads"
bash 00_env.sh > $L/00_env.log 2>&1
bash 01_region.sh > $L/01_region.log 2>&1
[[ -s $D/truth/truth_sv.vcf ]] || python3 02_somatic.py $D $CTG $SEED $N_SNV $N_IND $N_HP $N_DEL $N_INS $N_INV $N_DUP > $L/02_somatic.log 2>&1
samples | xargs -P $PAR -I{} bash -c "bash 03_simulate.sh {} > $L/03_simulate.{}.log 2>&1 || echo 'simulation {} FAILED'"
bash 04_normal.sh > $L/04_normal.log 2>&1
# raw tumor reads aligned once per sample: reference for the alternative reads before EC
samples | grep -v normal | xargs -P $PAR -I{} bash -c "source ./config.sh; source ./lib_sfs.sh; mkdir -p \$D/raw; [[ -s \$D/raw/{}.bam.bai ]] || align_tag \$D/reads/{}.fq.gz \$D/raw/{}.bam > $L/raw_align.{}.log 2>&1"
for s in $(samples | grep -v normal); do for e in $EC_SETS; do echo "$s $e"; done; done \
  | xargs -P $PAR -L1 bash -c "bash 05_cell.sh \$0 \$1 > $L/05_cell.\$0.\$1.log 2>&1 && python3 06_eval.py $D \$0 \$1 '$WEIGHTS' > $L/06_eval.\$0.\$1.log 2>&1 && python3 07_attrib.py $D \$0 \$1 > $L/07_attrib.\$0.\$1.log 2>&1 && python3 08_err.py $D \$0 \$1 $THREADS > $L/08_err.\$0.\$1.log 2>&1 || echo \"cell \$0 \$1 FAILED\""
mkdir -p $D/table
for t in small bp calls noise attr attr_background err; do
  f=$(ls $D/cells/*/*/eval/$t.tsv | head -1); head -1 $f > $D/table/$t.tsv
  for g in $D/cells/*/*/eval/$t.tsv; do tail -n +2 $g; done >> $D/table/$t.tsv
done
log "done: $D/table/ (small, bp, calls, noise, attr, attr_background, err)"
