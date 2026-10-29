#!/usr/bin/env bash
# Germline-removal cascade, extra row: index = GRCh38-GIABv3 + the curated HG008-N v6.3 assembly (one file: 46 chromosome
# haplotypes + 18 unassigned contigs). No normal reads in the index. It is the curated-assembly counterpart of the full
# germline set (reference + unitigs + missing reads, row 3 onwards): the claim to test is that the full set does better.
# Why v6.3 and not the FTP verkko 2.2: v6.3 is the same verkko 2.2 run after manual
# graph curation and polishing, 44/46 T2T scaffolds and 49 gaps; tumor and normal reads from the windows where verkko 2.2 is
# incomplete (acrocentric arms chr15/21/22, pericentromeric chr4/chr20, chrX:67-72 Mb) align to it at control divergence
# (refs/HG008N_v6.3/README.txt).
# CAVEATS for the text: (1) the v6.3 curation chose the paths through the normal graph that match the tumor assembly
# (except at somatic SVs) and scaffolded the acrocentric arms from the tumor, which can only favor this row; (2) the V0.5
# benchmark regions derive from verkko 2.2 and the V0.5 VNTR SVs from tumor-vs-v6.3 alignments; robustness claim only,
# never a cross-caller number.
# The uncapped SFS file is kept (search/tumor.ref_v63.sfs) to place the SFS longer than 500 bp.
# Needs: reads/tumor.ec.fa.gz from 00_reads.sh.
# Output: index/ref_v63.fmd, search/tumor.ref_v63.{sfs,dedup.cap500.sfs}, rows/R2w_refv63_{nofilters,normal}/, bench/summary.tsv
# Usage: nohup bash 06_v63.sh > logs/06_v63.log 2>&1 &
set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../../config.sh"
export PATH=$SV_BIN:$PATH
R=$HG008/reads_ec_align
D=$R/germline_cascade
GIABV3=$REFS/GRCh38_GIABv3.fasta.gz
V63=$HG008/refs/HG008N_v6.3/HG008N_v6.3.fasta.gz
REF=$R/ref/GRCh38_GIABv3.chr1-22X.fa
BAM=$R/align/tumor.ec.hp.bam
EXCL=$REFS/grch38_centromers.bed
THREADS=48; T=24
HPC="--hpc --hpc-cap 2"
NAME=R2w_refv63_nofilters; NAMEN=R2wb_refv63_normal
NBAM=$R/align/normal.ec.hp.bam
mkdir -p $D/index $D/search $D/time $D/logs $D/rows/$NAME $D/rows/$NAMEN
timed() { local name=$1; shift; /usr/bin/time -v -o $D/time/$name.txt "$@"; }
check_index()  { grep -qa "All done!" $1 && [[ -s $2 ]]        || { echo "index failed, see $1" >&2; exit 1; }; }
check_search() { grep -qa "All done!" $1 && grep -qv '^#' $2 || { echo "search failed, see $1" >&2; exit 1; }; }
[[ -s $D/reads/tumor.ec.fa.gz ]] || { echo "missing $D/reads/tumor.ec.fa.gz (00_reads.sh)" >&2; exit 1; }
gzip -t $V63 || { echo "corrupt $V63" >&2; exit 1; }

# 1. index: reference + v6.3, HPC cap 2
fmd=$D/index/ref_v63.fmd
if [[ ! -s $fmd ]]; then
  timed vespa_index_ref_v63 $VESPA index -d -m 4G -t $THREADS -p 0 $HPC $GIABV3 $V63 -o $fmd.part > $D/logs/vespa_index_ref_v63.log 2>&1
  check_index $D/logs/vespa_index_ref_v63.log $fmd.part
  mv $fmd.part $fmd; mv $fmd.part.hpc $fmd.hpc
fi

# 2. tumor search, then duplicates and SFS longer than 500 bp removed; the uncapped file is kept
raw=$D/search/tumor.ref_v63.sfs; sfs=$D/search/tumor.ref_v63.dedup.cap500.sfs
if [[ ! -s $raw ]]; then
  timed vespa_search_tumor_ref_v63 $VESPA search --index $fmd --fastx $D/reads/tumor.ec.fa.gz --bsize 10000 --threads $THREADS $HPC \
    > $raw.part 2> $D/search/vespa_search_tumor_ref_v63.log
  check_search $D/search/vespa_search_tumor_ref_v63.log $raw.part
  mv $raw.part $raw
fi
if [[ ! -s $sfs ]]; then
  python $LIB/dedup_sfs.py $raw $sfs.part --max-length 500 2> $D/search/dedup_sfs_ref_v63.log
  mv $sfs.part $sfs
fi

# 3. call without the normal, provenance, as the other rows (03_call_tables.sh)
d=$D/rows/$NAME
if [[ ! -s $d/calls.vcf ]]; then
  timed call_$NAME $VESPA call --reference $REF --bam $BAM --sfs $sfs --bed-exclusion $EXCL \
    --clip-clusters $d/calls.clipped_clusters --call-provenance $d/provenance.tsv \
    --min-cluster-weight 5 --min-sv-length 50 --min-mapq 60 --clipped --no-clipped-fallback --threads $T \
    > $d/calls.unsorted.vcf 2> $d/call.log
  bcftools sort $d/calls.unsorted.vcf -o $d/calls.vcf 2>/dev/null && rm -f $d/calls.unsorted.vcf
fi

# 3b. call with the normal (in-call filter, post-call flags)
dn=$D/rows/$NAMEN
if [[ ! -s $dn/calls.vcf ]]; then
  timed call_$NAMEN $VESPA call --reference $REF --bam $BAM --sfs $sfs --normal-contigs-bam $NBAM --bed-exclusion $EXCL \
    --clip-clusters $dn/calls.clipped_clusters --call-provenance $dn/provenance.tsv \
    --min-cluster-weight 5 --min-sv-length 50 --min-mapq 60 --clipped --no-clipped-fallback --threads $T \
    > $dn/calls.unsorted.vcf 2> $dn/call.log
  bcftools sort $dn/calls.unsorted.vcf -o $dn/calls.vcf 2>/dev/null && rm -f $dn/calls.unsorted.vcf
fi

# 4. breakpoint coverage
if [[ ! -s $d/sfs_clean.bed.gz ]]; then
  python $LIB/sfs_to_genome.py $sfs $BAM $d/sfs_supp.bed.gz $T 2> $d/sfs_to_genome.log
  zcat $d/sfs_supp.bed.gz | sort -k4,4 -k5,5n -k6,6n -S8G | awk -F'\t' '
    function flush(){ if(n==0) return; for(i=1;i<=n;i++) if(!hasexact || ex[i]==1) print buf[i] }
    { key=$4"\t"$5"\t"$6
      if(key!=prev){ flush(); n=0; hasexact=0; prev=key }
      n++; buf[n]=$0; ex[n]=$7; if($7==1) hasexact=1 }
    END{ flush() }' | sort -k1,1 -k2,2n -S8G | gzip > $d/sfs_clean.bed.gz.part
  mv $d/sfs_clean.bed.gz.part $d/sfs_clean.bed.gz; rm -f $d/sfs_supp.bed.gz
fi
FAI=$REF.fai TOL=50 python $REPRO/bp_coverage/bp_point_table_prov.py $d $d/provenance.tsv $NAME > $d/bp_point_prov_tol50.tsv

# 5. Minda, all rows
grep -q "^$NAME " $D/rows.tsv || echo "$NAME    $d/calls.vcf    0" >> $D/rows.tsv
grep -q "^R2wb_refv63_incall " $D/rows.tsv || echo "R2wb_refv63_incall    $dn/calls.vcf    0" >> $D/rows.tsv
grep -q "^R2wb_refv63_pass " $D/rows.tsv || echo "R2wb_refv63_pass    $dn/calls.vcf    1" >> $D/rows.tsv
bash $REPRO/germline_cascade/bench.sh
for f in $D/search/dedup_sfs_ref_utg.log $D/search/dedup_sfs_ref_v63.log $R/logs/dedup_sfs.log; do
  echo "$f: $(tr '\n' ' ' < $f)"; done
cat $d/bp_point_prov_tol50.tsv
echo done
