import sys
from collections import defaultdict

id, fai, hapblock = sys.argv[1:]

# build genome txt
L = list(range(1,23)) + ["X", "Y"]
if not "hs37d5" in fai:
    L = ["chr" + str(idx) for idx in L]

g = open("karyotype." + id + ".genome.txt", 'w')
g.write("chr\tstart\tend\n")
f = open(fai)
for line in f:
    chr, st, *_ = line.rstrip().split()
    if not chr in L: continue
    g.write(chr + "\t1\t" + st + "\n")
f.close()
g.close()

# build band txt

block = defaultdict(list)
chr, start, end = 0, 0, 0
f = open(hapblock)
for line in f:
    if line.startswith("*"): continue
    a = line.rstrip().split()

    if line.startswith("BLOCK"):
        if chr:
            block[chr].append(":".join([start,end]))
        line = next(f).strip()
        if line.startswith("*"): continue
        b = line.split()
        chr, start, end = b[3], b[4], b[4]

    else:
        end = a[4]
f.close()
for i,j in block.items():
    print(i,j)
sys.exit()

g = open("karyotype." + id + ".band.txt", 'w')
g.write("chr\tstart\tend\tname\tgieStain\n")
f = open("karyotype." + id + ".genome.txt")
f.readline()
for line in f:
    chr, st, ed = line.rstrip().split()
    pos, bn, bd = 1, 1, 1
    for hb in block[chr]:
        s, e = hb.split(":")
        s, e = int(s), int(e)
        if pos < s:
            g.write("\t".join([chr, str(pos), str(s), chr + "_" + str(bn), "gneg"]) + "\n")
            bn += 1
            pos = e
        elif pos == s:
            pos = e
        elif pos > e:
            continue
        col = "stalk" if bd % 2 == 0 else "gpos25"
        g.write("\t".join([chr, str(s), str(e), "band_" + str(bd), str(col)]) + "\n")
        bd += 1
    if pos < int(st):
        # $a[0]\t$pos\t$a[1]\t$a[0]\_$bn\tgneg\n
        g.write("\t".join([chr, str(pos), str(st), chr + "_" + str(bn), "gneg"]) + "\n")
        bn += 1

f.close()
g.close()