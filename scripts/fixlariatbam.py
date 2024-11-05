#!/usr/bin/python
## filter out reads != 100bp; add barcode to read id

import re,pickle,pysam,sys
from datetime import datetime

print(datetime.now())
fin, id, mapfile = sys.argv[1:]
fout = id + ".fixedlariat.bam"

D = {}
f = open(mapfile)
for line in f:
        a, b = line.split() # number, stlfr barcode
        D[a] = b
f.close()

print(datetime.now())

fin_h = pysam.AlignmentFile(fin, "rb")
fout_h = pysam.AlignmentFile(fout,'wb',header=fin_h.header)

for read in fin_h:
        # if len(read.seq) != 100:
        #         continue
        id = read.query_name    #ST-E0:0:SIMULATE:8:0:0:28782982
        # tag = read.get_tag("BX").replace("-1","")
        # new_id = id + "#" + D[tag] 
        a = id.split(":")[-1]
        new_id = "read" + a + "#" + D[a]
        read.query_name = new_id
        read.set_tag("BX",None)
        read.set_tag("AC",None)
        read.set_tag("XC",None)
        fout_h.write(read)
fin_h.close()
fout_h.close()
