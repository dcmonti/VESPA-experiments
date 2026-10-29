#!/bin/bash
# svviz2 2.0a3 (env svviz2) on the COLO829 truth set, tumor raw reads aligned with -s50; local realignment only
# (no bwa index of the reference), as the HG008 pilot (hg008/.../reads_ec_align/svviz2_pilot), report only.
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../../config.sh"
source $CONDA_SH && conda activate svviz2
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
cd $COLO829/reads_ec_align_s50/svviz2_vaf
/usr/bin/time -v svviz2 --ref $COLO829/reads_ec_align_s50/ref/GRCh38_GIABv3.chr1-22XY.fa --variants $COLO829/reads_ec_align_s50/svviz2_vaf/colo829_gt_svviz2.vcf --outdir $COLO829/reads_ec_align_s50/svviz2_vaf/reports \
    --only-realign-locally --no-render --no-dotplots "$COLO829/reads_ec_align_s50/align/tumor.raw.hp.bam,sequencer=pacbio"
