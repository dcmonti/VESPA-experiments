import sys
src, dst = sys.argv[1], sys.argv[2]
seen=set(); written=set(); cur=""; n_in=n_out=0
with open(src) as inf, open(dst,"w") as outf:
    for line in inf:
        s=line.rstrip("\n")
        if not s: continue
        if s.startswith("#"): outf.write(s+"\n"); continue
        n_in+=1; cols=s.split("\t")
        if cols[0]!="*": cur=cols[0]
        key=(cur,cols[1],cols[2])
        if key in seen: continue
        seen.add(key)
        if cur not in written: cols[0]=cur; written.add(cur)
        else: cols[0]="*"
        outf.write("\t".join(cols)+"\n"); n_out+=1
print(f"Input SFS entries: {n_in}")
print(f"Output SFS entries (dedup): {n_out}")
