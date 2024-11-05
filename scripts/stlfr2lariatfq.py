##barcode split stlfr reads -> lariat reads with stlfr barcode

import gzip, sys

barcodelist, fq1, fq2, outfq = sys.argv[1:]

h1 = gzip.open(fq1, 'rt')
h2 = gzip.open(fq2, 'rt')
wh = gzip.open(outfq,'w')

q = "F"*30
B = {}
f = open(barcodelist)
for line in f:
    id,seq = line.rstrip().split(",")
    B[id] = seq
f.close()

info = [] #id, bcseq, r1, r2, rq1, rq2
l = 0
for i,j in zip(h1, h2):
	l += 1
	if l % 4 == 1:
		if "0_0_0" in i:
			continue
		id = i.split()[0][:-2]
		if id.endswith("/1"): id = id.replace("/1","")
		bc = id.split("#")[1]
		a, b, c = bc.split("_")
		seq1, seq2, seq3 = B[a], B[b], B[c]
		seq = seq1 + seq2 + seq3 + "-1"
		info.append(id)
		info.append(seq)
	elif l % 4 == 2:
		if not info:
			continue
		info.append(i.strip())
		info.append(j.strip())
	elif l % 4 == 0:
		if len(info) != 4:
			continue
		info.append(i.strip())
		info.append(j.strip())
	if len(info) == 6:
		# print(info)
		wh.write(('\n'.join([info[0]+' 1:N:0:NAAGTGCT:0',info[2],info[4],info[3],info[5],info[1],q,'TTCACGCG','AAFFFKKK'])+'\n').encode())
		info = []
h1.close()
h2.close()
wh.close()
	
