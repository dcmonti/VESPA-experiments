#!/bin/bash
# svviz2 2.0a3 (env svviz2) on the HG008 V0.5 events of the clonal and subclonal benchmark, tumor raw reads
# (reads_ec_align/align/tumor.raw.hp.bam); settings as the COLO829 run: local realignment only, report only.
# Usage: run.sh <k>   (chunk k of chunks/, output reports/, log logs/chunk_<k>.log)
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../../config.sh"
source $CONDA_SH && conda activate svviz2
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
k=$1; mkdir -p $HG008/reads_ec_align/svviz2_vaf/reports $HG008/reads_ec_align/svviz2_vaf/logs $HG008/reads_ec_align/svviz2_vaf/tmp_$k; cd $HG008/reads_ec_align/svviz2_vaf/tmp_$k
/usr/bin/time -v svviz2 --ref $HG008/reads_ec_align/ref/GRCh38_GIABv3.chr1-22X.fa --variants $HG008/reads_ec_align/svviz2_vaf/chunks/chunk_$k.vcf --outdir $HG008/reads_ec_align/svviz2_vaf/reports \
    --only-realign-locally --no-render --no-dotplots "$HG008/reads_ec_align/align/tumor.raw.hp.bam,sequencer=pacbio" > $HG008/reads_ec_align/svviz2_vaf/logs/chunk_$k.log 2>&1
echo "chunk $k exit $?" >> $HG008/reads_ec_align/svviz2_vaf/logs/chunk_$k.log
