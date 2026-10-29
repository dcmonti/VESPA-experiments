#!/usr/bin/env bash
# Normal: full EC (hifiasm --het-cov NORMAL_HETCOV, as in production; the r_utg graph gives the
# contigs), alignment of the corrected reads + haplotag, index as in production:
# base = reference + contigs (cap 2), search of the normal against the base, reads with SFS (missing), full index.
set -euo pipefail
source "$(dirname "$0")/config.sh"; source $SRC/lib_sfs.sh
N=$D/normal; mkdir -p $N/ec $N/index
if [[ ! -s $N/ec/normal.ec.fa.gz ]]; then
  log "normal: hifiasm --het-cov $NORMAL_HETCOV"
  rm -rf $N/ec/work; mkdir -p $N/ec/work
  /usr/bin/time -v $HIFIASM -t $THREADS --write-ec --het-cov $NORMAL_HETCOV -o $N/ec/work/ec $D/reads/normal.fq.gz > $N/ec/hifiasm.log 2>&1
  gzip -1 -c $N/ec/work/ec.ec.fa > $N/ec/normal.ec.fa.gz
  awk '/^S/{print ">normal_"$2; print $3}' $(ls $N/ec/work/*r_utg.gfa | head -1) > $N/ec/contigs.fa
  rm -rf $N/ec/work
fi
[[ -s $N/normal.ec.bam.bai ]] || { log "normal: alignment"; align_tag $N/ec/normal.ec.fa.gz $N/normal.ec.bam; }
B=$N/index/base.fmd F=$N/index/full.fmd
[[ -s $B ]] || { log "base index"; $VESPA index -t $THREADS -d --hpc --hpc-cap 2 $D/ref/region.fa $N/ec/contigs.fa -o $B > $N/index/base.log 2>&1; }
if [[ ! -s $N/index/missing.fa ]]; then
  log "search of the normal against the base index"
  $VESPA search --index $B --fastx $N/ec/normal.ec.fa.gz --bsize 10000 --threads $THREADS --hpc --hpc-cap 2 > $N/index/normal.sfs 2> $N/index/normal.search.log
  python3 - $N/index/normal.sfs $N/ec/normal.ec.fa.gz $N/index/missing.fa <<'PY'
import sys, gzip
sfs, reads, out = sys.argv[1:4]
want=set(); last=""
for line in open(sfs):
    if line.startswith("#") or not line.strip(): continue
    r=line.split("\t")[0]
    if r=="*": r=last
    else: last=r
    want.add(r)
keep=False; n=0
with gzip.open(reads,"rt") as fh, open(out,"w") as fo:
    for line in fh:
        if line.startswith(">"):
            keep = line[1:].split()[0] in want; n += keep
        if keep: fo.write(line)
print(f"missing reads: {n}", file=sys.stderr)
PY
fi
[[ -s $F ]] || { log "full index"; $VESPA index -t $THREADS -d -i $B --hpc --hpc-cap 2 $N/index/missing.fa -o $F > $N/index/full.log 2>&1; }
log "normal: contigs $(grep -c '>' $N/ec/contigs.fa), missing $(grep -c '>' $N/index/missing.fa), index $(du -h $F | cut -f1)"
