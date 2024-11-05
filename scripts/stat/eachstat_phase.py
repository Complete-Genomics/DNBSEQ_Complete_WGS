## reseq:
# switch rate:          0.009980729514552672
# mismatch rate:        0.004383641570126749
# flat rate:            0.4541022747337955
# missing rate:         0.0006104136569882179
# phased count:         147351
# AN50:                 131078941.95383072
# N50:                  131495029
# max block snp frac:   0.6945950163160358

##new
# switch rate:        0.008333333333333333
# mismatch rate:      0.000630119722747322
# flat rate:          0.0008401596303297626
# phased count:       6566
# AN50:               0
# N50:                294
# num snps max blk:   79

import sys,gzip
from collections import OrderedDict

id, fin = sys.argv[1:]
print("chr\tswitch rate\tmismatch rate\tflat rate\tmissing rate\tphased count\tAN50\tN50\tmax block snp frac\tphasing rate")

num = 4
D = OrderedDict()
chrs = []
f = open(fin)
for line in f:
    line = line.strip()
    if not line: continue
    if "WARNING" in line: continue
    if not ":" in line:
        chr = line if not line.startswith("combine") else "chrAll"
        chrs.append(chr)
        D[chr] = []
    else:
        val = float(line.split()[-1])
        
        if val < 1:
            val = int((val * 10 ** num + 0.5)) / 10 ** num      #rounding to a fixed number of decimal places
            val = round(val, 4)
        else:
            val = int(val)
        D[chr].append(val)

f.close()

cnt, allratio = 0, 0
for chr in chrs:
    phased, het = 0, 0
    if chr == "chrAll": continue
    #hg002_prime.lariat.dv.chr3.hapblock.phased.VCF
    # print(id, chr)
    try:
        f = gzip.open(id + ".bwa.gatk." + chr + ".hapblock.phased.VCF.gz", 'rt')
    except:
        f = gzip.open(id + ".bwa.dv." + chr + ".hapblock.phased.VCF.gz", 'rt')
    for line in f:
        if line.startswith("#"): continue
        gt = line.rstrip().split()[9].split(":")[0]
        if "|" in gt:
            phased += 1
            het += 1
        if gt == "1/0" or gt == "0/1":
            het += 1
    ratio = round(phased/het, 4)
    D[chr].append(ratio)
    cnt += 1
    allratio += ratio
    f.close()
D["chrAll"].append(round(allratio/cnt,4))
if len(D["chrAll"]) == 8: D["chrAll"].insert(3, '0') # no missing rate
for chr, L in D.items():
    L = [str(i) for i in L]
    # print(chr, L, len(L))
    print(chr + "\t" + "\t".join(L))
