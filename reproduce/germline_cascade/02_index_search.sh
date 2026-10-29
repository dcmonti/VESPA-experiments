#!/usr/bin/env bash
# Germline-removal cascade, step 2: the two reduced indexes and the tumor search against each.
#   ref        GRCh38-GIABv3 only (all 218 contigs)                       -> rows 1 and 1b
#   ref_utg    GRCh38-GIABv3 + normal unitigs (the base index of 03)      -> row 2
# Same options as reproduce/hg008/03_vespa_index_search.sh: HPC cap 2, FMD format, then duplicate SFS and
# SFS longer than 500 bp removed. The full index (row 3 onwards) is the paper's index/full.fmd.
# PILOT=chrN restricts the tumor reads to those with primary alignment on chrN (to size the SFS set of the
# reference-only index before the genome-wide search); output then goes to search_pilot_chrN/.
# Output: index/{ref,ref_utg}.fmd, search/tumor.{ref,ref_utg}.dedup.cap500.sfs
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"
R=$HG008/reads_ec_align
D=$R/germline_cascade
GIABV3=$REFS/GRCh38_GIABv3.fasta.gz
THREADS=48
HPC="--hpc --hpc-cap 2"
INDEXES=${INDEXES:-"ref ref_utg"}
mkdir -p $D/index $D/time $D/logs
timed() { local name=$1; shift; /usr/bin/time -v -o $D/time/$name.txt "$@"; }
check_index()  { grep -qa "All done!" $1 && [[ -s $2 ]]        || { echo "index failed, see $1" >&2; exit 1; }; }
check_search() { grep -qa "All done!" $1 && grep -qv '^#' $2 || { echo "search failed, see $1" >&2; exit 1; }; }

TUMOR=$D/reads/tumor.ec.fa.gz; S=$D/search
if [[ -n ${PILOT:-} ]]; then
  S=$D/search_pilot_$PILOT; mkdir -p $S
  TUMOR=$S/tumor.$PILOT.ec.fa.gz
  [[ -s $TUMOR ]] || { $SAMTOOLS view -u -F 0x900 $R/align/tumor.ec.hp.bam $PILOT | $SAMTOOLS fasta - 2> $S/fasta.log | gzip -1 > $TUMOR.part; mv $TUMOR.part $TUMOR; }
fi
mkdir -p $S

for ix in $INDEXES; do
  fmd=$D/index/$ix.fmd
  if [[ ! -s $fmd ]]; then
    inputs=$GIABV3; [[ $ix == ref_utg ]] && inputs="$GIABV3 $D/unitigs/normal.r_utg.fa.gz"
    timed vespa_index_$ix $VESPA index -d -m 4G -t $THREADS -p 0 $HPC $inputs -o $fmd.part > $D/logs/vespa_index_$ix.log 2>&1
    check_index $D/logs/vespa_index_$ix.log $fmd.part
    mv $fmd.part $fmd; mv $fmd.part.hpc $fmd.hpc
  fi
  raw=$S/tumor.$ix.sfs; final=$S/tumor.$ix.dedup.cap500.sfs
  [[ -s $final ]] && continue
  timed vespa_search_tumor_$ix${PILOT:+_$PILOT} $VESPA search --index $fmd --fastx $TUMOR --bsize 10000 --threads $THREADS $HPC \
    > $raw.part 2> $S/vespa_search_tumor_$ix.log
  check_search $S/vespa_search_tumor_$ix.log $raw.part
  mv $raw.part $raw
  $PYTHON $LIB/dedup_sfs.py $raw $final.part --max-length 500 2> $S/dedup_sfs_$ix.log
  mv $final.part $final
  rm -f $raw
done
for ix in $INDEXES; do echo "$ix: $(grep -vc '^#' $S/tumor.$ix.dedup.cap500.sfs) SFS"; done
