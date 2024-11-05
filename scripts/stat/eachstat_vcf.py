import sys

sample = sys.argv[1]

col1 = ["Sample", "Total_SNP", "dbSNP_rate", "Novel_SNP", "Novel_SNP_Rate", "Ti/Tv", "Total_INDEL", "dbINDEL_Rate"]
col2 = [sample]

allsnp = dbsnp = novel = ti = tv = 0
f = open("snps.vcf")
for line in f:
    if line.startswith("#"): continue
    line = line.strip().split()
    if line[6] == "RefCall": continue
    allsnp += 1
    if line[2] == ".": 
        novel += 1
    else:
        dbsnp += 1
    ref, alt = line[3], line[4].split(",")
    if len(alt) > 1: continue
    alt = alt[0]
    if  (ref == 'A' and alt == 'G') or (ref == 'G' and alt == 'A') or \
        (ref == 'C' and alt == 'T') or (ref == 'T' and alt == 'C'):
        ti += 1
    else:
        tv += 1
f.close()
dbsnprate = int(dbsnp / allsnp * 10000 + 0.5 ) / 100 #int($stat[1] / $stat[0] * 10000 + 0.5) / 100
novelrate = int((100 - dbsnprate) * 100 + 0.5 ) / 100 #int((100 - $stat[1]) * 100 + 0.5) / 100
tt = int(ti/tv * 100 + 0.5) /100 #int($ti / $tv * 100 + 0.5) / 100

allsnp = str(allsnp)
novel = str(novel)
dbsnprate = str(dbsnprate) + "%"
novelrate = str(novelrate) + "%"
tt = str(tt)

col2 += [allsnp, dbsnprate, novel, novelrate, tt]

indel = dbindel = 0
f = open("indels.vcf")
for line in f:
    if line.startswith("#"): continue
    line = line.strip().split()
    if line[6] == "RefCall": continue
    indel += 1
    if line[2] != ".": 
        dbindel += 1
f.close()
dbindelrate = int(dbindel / indel * 10000 + 0.5 ) / 100
dbindelrate = str(dbindelrate) + "%"
indel = str(indel)
col2 += [indel, dbindelrate]
for i in range(len(col1)):
    print(col1[i] + "\t" + col2[i])