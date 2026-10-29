#!/usr/bin/env bash
# HG008: VESPA index and search, from the error-corrected reads and the normal unitigs written by 01_reads_ec_align.sh.
#
# The index holds the reference and the germline sequence of the normal, so that the tumor search
# reports only sequence absent from both (sample-specific strings, SFS):
#   1. base index: GRCh38-GIABv3 (all 218 contigs) + normal unitigs
#   2. search the normal reads against the base index
#   3. normal reads that still produce SFS: they carry germline sequence the unitigs miss
#   4. full index: base index + those reads
#   5. search the tumor reads against the full index, remove duplicate SFS and SFS longer than 500 bp
# Index and search homopolymer-compress reference and reads, with runs capped at 2.
#
# Output, in $OUT:
#   index/full.fmd, index/full.fmd.hpc   full index
#   search/tumor.dedup.cap500.sfs        SFS given to VESPA call
#   time/<step>.txt                      GNU time (/usr/bin/time -v) report of each step
# Everything else this script reads or writes is deleted at the end: base index, normal SFS, normal reads
# added to the index, raw tumor SFS, error-corrected reads, normal unitigs.
#
# Usage: bash 03_vespa_index_search.sh
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"

OUT=${OUT:-$HG008/reads_ec_align}
GIABV3=$REFS/GRCh38_GIABv3.fasta.gz

THREADS=48
HPC="--hpc --hpc-cap 2"
MAX_SFS_LENGTH=500

mkdir -p $OUT/index $OUT/search $OUT/logs $OUT/time

# timed NAME COMMAND...: run COMMAND under GNU time, report in $OUT/time/NAME.txt.
timed() { local name=$1; shift; /usr/bin/time -v -o $OUT/time/$name.txt "$@"; }

TUMOR_EC=$OUT/ec/tumor.ec.fa.gz
NORMAL_EC=$OUT/ec/normal.ec.fa.gz
UNITIGS=$OUT/ec/normal.r_utg.fa.gz
BASE=$OUT/index/base.fmd
FULL=$OUT/index/full.fmd
NORMAL_SFS=$OUT/search/normal.sfs
MISSING=$OUT/search/normal.missing.fa.gz
TUMOR_SFS=$OUT/search/tumor.sfs
FINAL_SFS=$OUT/search/tumor.dedup.cap$MAX_SFS_LENGTH.sfs

[[ -s $FINAL_SFS ]] && { echo "already done: $FINAL_SFS"; exit 0; }

# VESPA writes "All done!" at the end of an index or a search, but also when an input is missing:
# the output is checked as well (a non-empty index, at least one SFS).
check_index()  { grep -qa "All done!" $1 && [[ -s $2 ]]        || { echo "index failed, see $1" >&2; exit 1; }; }
check_search() { grep -qa "All done!" $1 && grep -qv '^#' $2 || { echo "search failed, see $1" >&2; exit 1; }; }

# 1. Base index. -d writes the FMD format, the one VESPA search reads.
if [[ ! -s $BASE.hpc ]]; then
  timed vespa_index_base $VESPA index -d -m 4G -t $THREADS -p 0 $HPC $GIABV3 $UNITIGS -o $BASE.part > $OUT/logs/vespa_index_base.log 2>&1
  check_index $OUT/logs/vespa_index_base.log $BASE.part
  mv $BASE.part $BASE
  mv $BASE.part.hpc $BASE.hpc
fi

# 2. Normal search against the base index.
if [[ ! -s $NORMAL_SFS ]]; then
  timed vespa_search_normal $VESPA search --index $BASE --fastx $NORMAL_EC --bsize 10000 --threads $THREADS $HPC \
    > $NORMAL_SFS.part 2> $OUT/logs/vespa_search_normal.log
  check_search $OUT/logs/vespa_search_normal.log $NORMAL_SFS.part
  mv $NORMAL_SFS.part $NORMAL_SFS
fi

# 3. Normal reads with at least one SFS.
if [[ ! -s $MISSING ]]; then
  timed extract_missing_reads $PYTHON $LIB/extract_missing_reads.py $NORMAL_SFS $NORMAL_EC $MISSING.part.gz 2> $OUT/logs/extract_missing_reads.log
  mv $MISSING.part.gz $MISSING
fi

# 4. Full index: the reads are appended to the base index (-i).
if [[ ! -s $FULL.hpc ]]; then
  timed vespa_index_full $VESPA index -d -m 4G -t $THREADS -p 0 $HPC -i $BASE $MISSING -o $FULL.part > $OUT/logs/vespa_index_full.log 2>&1
  check_index $OUT/logs/vespa_index_full.log $FULL.part
  mv $FULL.part $FULL
  mv $FULL.part.hpc $FULL.hpc
fi

# 5. Tumor search against the full index, then duplicates and long SFS removed.
if [[ ! -s $TUMOR_SFS ]]; then
  timed vespa_search_tumor $VESPA search --index $FULL --fastx $TUMOR_EC --bsize 10000 --threads $THREADS $HPC \
    > $TUMOR_SFS.part 2> $OUT/logs/vespa_search_tumor.log
  check_search $OUT/logs/vespa_search_tumor.log $TUMOR_SFS.part
  mv $TUMOR_SFS.part $TUMOR_SFS
fi
timed dedup_sfs $PYTHON $LIB/dedup_sfs.py $TUMOR_SFS $FINAL_SFS.part --max-length $MAX_SFS_LENGTH 2> $OUT/logs/dedup_sfs.log
mv $FINAL_SFS.part $FINAL_SFS

rm -f $BASE $BASE.hpc
rm -f $NORMAL_SFS $MISSING $TUMOR_SFS $TUMOR_EC $NORMAL_EC $UNITIGS
