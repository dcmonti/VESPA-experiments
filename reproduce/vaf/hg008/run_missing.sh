#!/bin/bash
# Final pass of the HG008 svviz2 run on the groups without a report (make_missing.py); settings as run.sh.
# Usage: run_missing.sh <k>   (chunk k of missing/, output reports/, log logs/missing_<k>.log)
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../../config.sh"
source $CONDA_SH && conda activate svviz2
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
k=$1; mkdir -p $HG008/reads_ec_align/svviz2_vaf/tmp_missing_$k; cd $HG008/reads_ec_align/svviz2_vaf/tmp_missing_$k
/usr/bin/time -v svviz2 --ref $HG008/reads_ec_align/ref/GRCh38_GIABv3.chr1-22X.fa --variants $HG008/reads_ec_align/svviz2_vaf/missing/chunk_$k.vcf --outdir $HG008/reads_ec_align/svviz2_vaf/reports \
    --only-realign-locally --no-render --no-dotplots "$HG008/reads_ec_align/align/tumor.raw.hp.bam,sequencer=pacbio" > $HG008/reads_ec_align/svviz2_vaf/logs/missing_$k.log 2>&1
echo "missing $k exit $?" >> $HG008/reads_ec_align/svviz2_vaf/logs/missing_$k.log
