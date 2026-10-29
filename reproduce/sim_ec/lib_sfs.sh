# Shared functions: alignment + haplotag, search + dedup + cut at 500 bp (as in run_hpc_caps.sh).
align_tag(){ # fastx out.bam  -> minimap2 -s100 on the simulated reference, HP from the read name
  local IN=$1 OUT=$2; local T=${OUT%.bam}
  minimap2 -ax map-hifi -s100 --secondary=no --eqx --cs --MD -t $THREADS $D/ref/region.fa $IN 2> $T.minimap2.log \
    | samtools sort -@ 4 -m 2G -o $T.untagged.bam -
  # name <sample>.<component>.<n>: hap2 -> H2, hap1 / hap1_T1 / hap1_T2 -> H1; UNTAG_FRAC untagged.
  # The component is the SECOND-TO-LAST field: the sample name can contain a dot (t60_p0.3);
  # taking the second field tagged all hap2 reads of the tumor as H1 (fixed on 2026-09-29).
  samtools view $T.untagged.bam | cut -f1 | sort -u | awk -v f=$UNTAG_FRAC -v s=$SEED -v c=$CTG 'BEGIN{srand(s); print "#readname\thaplotype\tphaseset\tchromosome"}
    { n=split($1,a,"."); h=(a[n-1]=="hap2")?"H2":"H1"; if(rand()<f) h="none"; print $1"\t"h"\t"(h=="none"?"none":1)"\t"c }' > $T.haplotags.tsv
  python $TAGBAM $T.untagged.bam $T.haplotags.tsv $OUT --threads 4 --with-ps > $T.tag.log 2>&1
  samtools index $OUT && rm -f $T.untagged.bam
}
search_sfs(){ # index fastx out_prefix -> <out>.dedup.cap500.sfs
  local IDX=$1 IN=$2 P=$3
  $VESPA search --index $IDX --fastx $IN --bsize 10000 --threads $THREADS --hpc --hpc-cap 2 > $P.sfs 2> $P.search.log
  grep -qa "All done!" $P.search.log || { echo "search not completed ($P.search.log)" >&2; return 1; }
  python3 $SRC/dedup_sfs.py $P.sfs $P.dedup.sfs > $P.dedup.log
  awk -v MAX=500 'BEGIN{FS=OFS="\t"} /^#/{print; next}
    { if($1!="*"){cur=$1; seen=0}
      if($3<=MAX){ if(!seen){$1=cur; seen=1} else {$1="*"} print } }' $P.dedup.sfs > $P.dedup.cap500.sfs
  rm -f $P.sfs $P.dedup.sfs
}
