import gzip, sys
from collections import defaultdict

fin = sys.argv[1]

D = defaultdict(int)
l = 0
all_reads, split_reads = 0, 0

f = gzip.open(fin, 'rt')
for line in f:
    l += 1
    if l % 4 != 1: continue
    all_reads += 1
    id = line.rstrip().split()[0].split("#")[1]
    D[id] += 1
    if not "0_0_0" in line: 
        split_reads += 1
f.close()

#print(round(split_reads/all_reads, 2))
for _ in range(4): print(".")
ratio = float(split_reads)/float(all_reads) * 100
line = f"Reads_pair_num(after split) = {split_reads}({ratio}%)"
print(line)
for bc, cnt in D.items():
    print("\t".join([".",str(cnt), str(bc)]))

