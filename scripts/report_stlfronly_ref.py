# id, aligner, varcaller, stlfrbam avg depth, 

import sys,os,re, subprocess


"""
[biancai@cngb-xcompute-0-17 G400_ECR6_stLFR-1]$ head 01.filter/G400_ECR6_stLFR-1/split_stat_read1.log
Barcode_types = 1536*1536*1536=3623878656
Barcode_types_with_mismatch = 62976*62976*62976=249761340850176(1)
Real_Barcode_types = 33574842
Reads_pair_num = 818923261
Reads_pair_num(after split) = 753672975(92.032186%)
0       65250286        0_0_0
8282618 1       1000_1000_1064
11119045        1       1000_1000_1096

"""
id, aligner, varcaller, varstats, het, splitLog, lfr, aligncatstlfr,phase, stlfrbamdepth= sys.argv[1:]
L = [id, aligner, varcaller, stlfrbamdepth] ##
##snp indel count
f = open(varstats)
for line in f:
	if not line.startswith("SN"): continue
	if "number of SNPs" in line: 
		snp = line.rstrip().split()[-1]
	elif "number of indels" in line:
		indel = line.rstrip().split()[-1]
		break
f.close()

f = open(het)
hetsnp, hetindel = f.readline().rstrip().split()
f.close()
# process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# stdout, stderr = process.communicate()
# counts = stdout.decode().strip().replace("\n","\t")
# hetsnp, hetindel = counts.split()

for i in [snp, hetsnp, indel, hetindel]: L.append(i) ##



##barcode split rate (pe reads num/pe reads num after split)
f = open(splitLog)
try:
	for _ in range(4): f.readline()
	bcsplitrate = re.split('\(|\%\)', f.readline().split()[-1])[1]
	bcsplitrate = round(float(bcsplitrate), 2)
except:
	bcsplitrate = f.readline().strip()
L.append(bcsplitrate)								##
f.close()



##lfr num and avg len
"""
good-lfr        2124431
Average-lfr-length      184375
Average-lfr-readpair    7
(valid_lfr/valid_barcode)       1.35092

"""
f = open(lfr)
for line in f:
	a = line.rstrip().split()[1]
	if "good" in line:
		lfrnum = a
	elif "length" in line:
		avglen = a
	elif "readpair" in line:
		avgfragreadcount = a
	else:
		lfrperbc = a
f.close()
L.append(lfrnum)						##
L.append(avglen)						##		

## map stats #sample, map rate, PE map rate, meanis, duprate, cov1pct cov10pct
# stLFR
# 96.27%
# 92.57%
# 300
# 0.13
# 0.971404 0.962186

f = open(aligncatstlfr)
for _ in range(2): f.readline()
pemapratestlfr = f.readline().rstrip()
for _ in range(2): f.readline()
cov10long = f.readline().rstrip().split()[-1]
f.close()

L.append(pemapratestlfr)				##
L.append(cov10long)						##

##phase
f = open(phase)
allhetsnp, phasedhetsnp, allhetindel, phasedhetindel, fbn50, fbnum = f.readline().rstrip().split()
f.close()


for i in [phasedhetsnp, phasedhetindel, fbnum, fbn50]:
	L.append(i)							##							

for i in L:
	print(i)

# for i in [id, snp, hetsnp, indel, hetindel, bcsplitrate, lfrnum, avglen, pemapratestlfr, pemapratepf, cov10long, cov10short, phasedhetsnp, phasedhetindel, fbnum, fbn50, genecov, cmrggenecov]:
