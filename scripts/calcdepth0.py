import sys
from collections import defaultdict

fin = sys.argv[1]

euchs = ["chr" + str(_) for _ in range(1,23)]

#poss = 0
depthD = defaultdict(int)
f = open(fin)
for line in f:
    line = line.rstrip().split()
    try:
        chr, pos, depth = line
    except:
        continue
    # only euchromatic
#    if not chr in euchs:
#        continue
    depth = int(depth)
    if not depth:
        continue    
#    poss += 1
    depthD[depth] += 1  #depth: pos count
f.close()

totalcovbase, covbase, covbase10, covbase20, covbase30 = 0, 0, 0, 0, 0
for depth in depthD:
    count = depthD[depth]
    covbase += count
    totalcovbase += depth*count
    if depth >= 10:
        covbase10 += count
    if depth >= 20:
        covbase20 += count
    if depth >= 30:
        covbase30 += count
    

print(covbase, covbase10)
