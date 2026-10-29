#!/bin/bash
# svviz2 input from the COLO829 truth set of the bench (bench/min50/gt.vcf), 2026-10-27:
#  - the 3 INS are dropped: symbolic <INS> without sequence, which svviz2 cannot build an alt allele for
#  - <DUP> -> <DUP:TANDEM>, the only duplication svviz2 accepts
#  - BND: the anchor base inside ALT is set to REF (13 records keep the GRCh37 base after the liftover);
#    svviz2 uses only position and orientation, and stops on ALT base != REF ("complex event")
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../../config.sh"
G=$COLO829/reads_ec_align_s50/bench/min50/gt.vcf
O=$COLO829/reads_ec_align_s50/svviz2_vaf
(grep '^#' $G
 grep -v '^#' $G | awk -F'\t' 'BEGIN{OFS="\t"}
   $8~/SVTYPE=INS/ {next}
   $8~/SVTYPE=DUP/ {$5="<DUP:TANDEM>"}
   $8~/SVTYPE=BND/ { if (match($5,/^[][][^][]*[][]/)) $5=substr($5,1,RLENGTH) $4; else if (match($5,/[][][^][]*[][]$/)) $5=$4 substr($5,RSTART) }
   {print}') > $O/colo829_gt_svviz2.vcf
