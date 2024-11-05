#!/usr/bin/env python
####usage: python3 stlfr2lariat_v2.py sample_S1_L001_R1_001.fastq.gz sample_S1_L001_R2_001.fastq.gz input_4_lariat.fq.gz

import sys
import gzip

infile1=gzip.open(sys.argv[1])
infile2=gzip.open(sys.argv[2])
outfile=gzip.open(sys.argv[3],'w')

dic={}
line_num=1

info=[]

for i,j in zip(infile1,infile2):
    i=i.decode().strip()
    j=j.decode().strip()
    if line_num % 4==1:
        id=i.split(' ')[0]
        info.append(id)
    elif line_num % 4==2:
        fa1=i[16:]
        fa2=j
        barcode=i[0:16]
        info.append(fa1);info.append(fa2);info.append(barcode)
    elif line_num % 4==0:
        fq1=i[16:]
        fq2=j
        barcodeq=i[0:16]
        info.append(fq1);info.append(fq2);info.append(barcodeq)
    else:
        pass
    line_num+=1
    if len(info)==7:
        outfile.write(('\n'.join([info[0]+' 1:N:0:NAAGTGCT:0',info[1],info[4],info[2],info[5],info[3]+'-1,'+info[3],info[6],'TTCACGCG','AAFFFKKK'])+'\n').encode())
        info=[]
    else:
        pass

outfile.close()
