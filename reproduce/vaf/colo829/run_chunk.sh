#!/bin/bash
# svviz2 on one chunk of the COLO829 events without a report (make_chunks.py); settings as run.sh.
# Usage: run_chunk.sh <k>   (output reports/, log logs/chunk_<k>.log)
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../../config.sh"
source $CONDA_SH && conda activate svviz2
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
k=$1; mkdir -p $COLO829/reads_ec_align_s50/svviz2_vaf/tmp_$k; cd $COLO829/reads_ec_align_s50/svviz2_vaf/tmp_$k
/usr/bin/time -v svviz2 --ref $COLO829/reads_ec_align_s50/ref/GRCh38_GIABv3.chr1-22XY.fa --variants $COLO829/reads_ec_align_s50/svviz2_vaf/chunks/chunk_$k.vcf --outdir $COLO829/reads_ec_align_s50/svviz2_vaf/reports \
    --only-realign-locally --no-render --no-dotplots "$COLO829/reads_ec_align_s50/align/tumor.raw.hp.bam,sequencer=pacbio" > $COLO829/reads_ec_align_s50/svviz2_vaf/logs/chunk_$k.log 2>&1
echo "chunk $k exit $?" >> $COLO829/reads_ec_align_s50/svviz2_vaf/logs/chunk_$k.log
