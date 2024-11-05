import sys
from collections import defaultdict

pfbed, mergebed, id = sys.argv[1:]

D = defaultdict(float)
f = open(pfbed)
for line in f:
    line = line.rstrip().split()
    if len(line) != 8: continue
    chr, st, ed, gene, depth, cnt, l, ratio = line
    if chr == "all" or int(depth) < 10 : continue
    D[gene] += float(ratio)
f.close()

g = open(id + ".cmrg.pf.cov", 'w')
for gene, cov in D.items():
    cov = str(cov)
    g.write(gene + "\t" + cov + "\n")
g.close()

D = defaultdict(float)
f = open(mergebed)
for line in f:
    line = line.rstrip().split()
    if len(line) != 8: continue
    chr, st, ed, gene, depth, cnt, l, ratio = line
    if chr == "all" or int(depth) < 10 : continue
    D[gene] += float(ratio)
f.close()

g = open(id + ".cmrg.merge.cov", 'w')
for gene, cov in D.items():
    cov = str(cov)
    g.write(gene + "\t" + cov + "\n")
g.close()