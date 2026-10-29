#!/usr/bin/env bash
# Fix of 2026-09-29: the HP tags of the tumor were wrong (all hap2 reads with HP=1, see
# lib_sfs.sh). For one cell: rebuilds the tag list from the read names (second-to-last field), retags the
# corrected, already aligned BAM (tag_bam_from_list.py removes the existing tags), deletes dump, calls and evaluation
# and reruns them (05_cell.sh skips EC, alignment and search, which do not depend on the tags), then 06_eval.py and
# 07_attrib.py. USAGE: bash 08_retag_recall.sh t60_p0.3 fnr0.10
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"; source ./config.sh
S=${1:?sample}; EC=${2:?setting}; W=$D/cells/$S/$EC; B=$W/tumor.ec.bam; T=${B%.bam}
samtools view $B | cut -f1 | sort -u | awk -v f=$UNTAG_FRAC -v s=$SEED -v c=$CTG 'BEGIN{srand(s); print "#readname\thaplotype\tphaseset\tchromosome"}
  { n=split($1,a,"."); h=(a[n-1]=="hap2")?"H2":"H1"; if(rand()<f) h="none"; print $1"\t"h"\t"(h=="none"?"none":1)"\t"c }' > $T.haplotags.tsv
python $TAGBAM $B $T.haplotags.tsv $T.retag.bam --threads 4 --with-ps > $T.retag.log 2>&1
mv $T.retag.bam $B && samtools index $B
rm -rf $W/prov $W/prov_w1 $W/call_w* $W/eval
bash 05_cell.sh $S $EC
python3 06_eval.py $D $S $EC "$WEIGHTS"
python3 07_attrib.py $D $S $EC
log "$S $EC: retagged and rerun ($(samtools view -F 0x900 $B | awk '{for(i=12;i<=NF;i++) if($i~/^HP:i:/) c[substr($i,6)]++} END{printf "HP1 %d, HP2 %d", c[1], c[2]}'))"
