#!/usr/bin/env bash
# Germline-removal cascade, step 0: the error-corrected reads, recovered from the EC alignments.
# 03_vespa_index_search.sh deleted ec/{tumor,normal}.ec.fa.gz after the search. The EC BAMs hold every
# corrected read once as a primary or unmapped record (-F 0x900); samtools fasta writes reverse-strand
# records back in their original orientation, so the sequences are the ones the search read.
# Output: reads/tumor.ec.fa.gz, reads/normal.ec.fa.gz (the normal is used only to check 01_unitigs.sh).
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"
R=$HG008/reads_ec_align
D=$R/germline_cascade
mkdir -p $D/reads $D/time
for s in tumor normal; do
  out=$D/reads/$s.ec.fa.gz
  [[ -s $out ]] && continue
  /usr/bin/time -v -o $D/time/fasta.$s.txt bash -o pipefail -c \
    "$SAMTOOLS fasta -@ 8 -F 0x900 $R/align/$s.ec.hp.bam 2> $D/reads/$s.fasta.log | gzip -1 > $out.part"
  mv $out.part $out
done
