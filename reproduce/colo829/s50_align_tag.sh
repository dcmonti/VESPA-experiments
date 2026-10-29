#!/usr/bin/env bash
# COLO829, second alignment: the same reads aligned with minimap2 -ax map-hifi -s50, and the haplotype tags of the
# first run transferred by read name.
#
# Why: the map-hifi preset sets the minimum chaining score to 200 (-s200, match score 1), so a terminal segment of a
# read shorter than about 200 bp gets no supplementary alignment. On COLO829, 9 of the 62 truth events (mostly
# translocations of the chr3/chr10/chr12 and chr15/chr20 clusters, and the 289 bp inversion on chr6) are found by a
# published callset and by none of the five callers on the default alignment; at the chr3-chr12 junction the PacBio
# BAM (pbmm2) has 63 split reads and the default alignment none. -s50 is the setting of Qin & Li on these reads, and
# keeps supplementary segments down to about 50 bp, as pbmm2 does (--min-length 50). On HG008 the default alignment
# loses 1 event of 172, so HG008 keeps the default alignment only.
#
# Input (from the first run, $OUT_DEFAULT, written by 01_reads_ec_align.sh and 02_haplotag.sh):
#   align/{tumor,normal}.{raw,ec}.hp.bam   the reads are taken back from these BAMs (primary and unmapped records,
#                                          -F 0x900), so both alignments see exactly the same reads; the
#                                          error-corrected reads, deleted after the first run, are recovered this way
#   haplotag/{tumor,normal}.haplotags.tsv  whatshap read lists, applied by read name (no new phasing or haplotagging:
#                                          the read set and the phased VCF of the normal are the same)
#   ref/, haplotag/normal.phased.vcf.gz, search/tumor.dedup.cap500.sfs   linked into $OUT, so that 04_call.sh and
#                                          05_bench.sh run unchanged with OUT=$OUT (the SFS do not depend on the
#                                          alignment: the search reads the reads, not the BAM)
#
# Output, in $OUT:
#   align/{tumor,normal}.{raw,ec}.hp.bam   tagged alignments, indexed
#   time/<step>.txt                        GNU time (/usr/bin/time -v) report of each step
#
# Usage: bash s50_align_tag.sh            (s50_run_all.sh runs it, then 04_call.sh and 05_bench.sh)
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"

OUT_DEFAULT=${OUT_DEFAULT:-$COLO829/reads_ec_align}
OUT=${OUT:-$COLO829/reads_ec_align_s50}

TAG_BAM=$LIB/tag_bam_from_list.py

THREADS=48
PARALLEL=${PARALLEL:-2}   # alignments at a time (tumor and normal of the same reads); ~67 GB of memory each

mkdir -p $OUT/align $OUT/logs $OUT/time $OUT/search
REF=$OUT_DEFAULT/ref/GRCh38_GIABv3.chr1-22XY.fa
[[ -e $OUT/ref ]]                        || ln -s $OUT_DEFAULT/ref $OUT/ref
[[ -e $OUT/haplotag ]]                   || ln -s $OUT_DEFAULT/haplotag $OUT/haplotag
[[ -e $OUT/search/tumor.dedup.cap500.sfs ]] || ln -s $OUT_DEFAULT/search/tumor.dedup.cap500.sfs $OUT/search/tumor.dedup.cap500.sfs

# timed NAME COMMAND...: run COMMAND under GNU time, report in $OUT/time/NAME.txt.
timed() { local name=$1; shift; /usr/bin/time -v -o $OUT/time/$name.txt "$@"; }

# one alignment: reads of the first-run BAM -> minimap2 -s50 -> sorted BAM -> HP/PS tags -> index
align_tag() { # sample reads
  local sample=$1 reads=$2
  local src=$OUT_DEFAULT/align/$sample.$reads.hp.bam
  local bam=$OUT/align/$sample.$reads.bam hp=$OUT/align/$sample.$reads.hp.bam
  [[ -s $hp.bai ]] && return 0
  # raw reads have qualities (FASTQ), error-corrected reads are FASTA, as in the first run
  local fmt=fastq; [[ $reads == ec ]] && fmt=fasta
  if [[ ! -s $bam ]]; then
    timed align.$sample.$reads bash -o pipefail -c "$SAMTOOLS $fmt -@ 8 -F 0x900 $src 2> $OUT/logs/$fmt.$sample.$reads.log \
      | $MINIMAP2 -ax map-hifi -s50 -t $THREADS $REF - 2> $OUT/logs/minimap2.$sample.$reads.log \
      | $SAMTOOLS sort -@ 16 -m 4G -T $OUT/align/$sample.$reads.sort -o $bam.part -"
    mv $bam.part $bam
  fi
  timed tag.$sample.$reads $PYTHON $TAG_BAM $bam $OUT_DEFAULT/haplotag/$sample.haplotags.tsv $hp.part --threads 16 --with-ps \
    2> $OUT/logs/tag.$sample.$reads.log
  mv $hp.part $hp
  timed index.$sample.$reads.hp $SAMTOOLS index -@ 16 $hp
  rm -f $bam
}

for reads in raw ec; do
  pids=()
  for sample in tumor normal; do
    align_tag $sample $reads > $OUT/logs/align_tag.$sample.$reads.log 2>&1 &
    pids+=($!)
    [[ $PARALLEL -ge 2 ]] || wait ${pids[-1]}
  done
  for p in "${pids[@]}"; do wait $p || { echo "alignment $reads failed, see $OUT/logs/" >&2; exit 1; }; done
done

# read counts: the second alignment must hold the same reads as the first (primary + unmapped records)
for sample in tumor normal; do for reads in raw ec; do
  a=$($SAMTOOLS view -@ 8 -c -F 0x900 $OUT_DEFAULT/align/$sample.$reads.hp.bam)
  b=$($SAMTOOLS view -@ 8 -c -F 0x900 $OUT/align/$sample.$reads.hp.bam)
  s=$($SAMTOOLS view -@ 8 -c -f 0x800 $OUT/align/$sample.$reads.hp.bam)
  echo -e "$sample\t$reads\treads_default=$a\treads_s50=$b\tsupplementary_s50=$s"
  [[ $a == $b ]] || { echo "read count differs for $sample $reads" >&2; exit 1; }
done; done | tee $OUT/align/read_counts.tsv
