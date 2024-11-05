import sys,os
import re

fflagstat, fstats, fdepth, finsertsize, sample = sys.argv[1:]
L = []
"""flagstat:
1737882 + 0 in total (QC-passed reads + QC-failed reads)
204 + 0 secondary
4996 + 0 supplementary
479 + 0 duplicates
1737158 + 0 mapped (99.96% : N/A)
1732682 + 0 paired in sequencing
866806 + 0 read1
865876 + 0 read2
1700317 + 0 properly paired (98.13% : N/A)
1726960 + 0 with itself and mate mapped
4998 + 0 singletons (0.29% : N/A)
15892 + 0 with mate mapped to a different chr
10125 + 0 with mate mapped to a different chr (mapQ>=5)
"""
# f = open(fstats)
# for line in f:
#     if not line.startswith("SN"):
#         continue
#     line =line.strip()
#     if "1st fragments" in line:
#         reads = int(line.split()[0]) * 2
#     if "reads duplicated" in line:
#         dupreads = int(line.split()[3])
# f.close()

# f = open(fflagstat)
dupreads = None
maprate = None
pemapreads = None
pemaprate = None

with open(fflagstat, 'r') as file:
    for line in file:
        match = re.search(r'(\d+) \+ \d+ duplicates', line)
        if match:
            dupreads = int(match.group(1))
        else:
            match = re.search(r'\+ \d+ mapped \(([\d.]+%) : \S+\)', line)
            if match:
                maprate = match.group(1)
                L.append(maprate)
            else:
                match = re.search(r'(\d+) \+ \d+ properly paired \(([\d.]+%) : \S+\)', line)
                if match:
                    pemapreads = int(match.group(1))
                    pemaprate = match.group(2)
                    L.append(pemaprate)

# for _ in range(3): f.readline()
# dupreads = int(f.readline().split()[0])
# maprate = f.readline().split()[4].replace("(","")
# L.append(maprate)
# for _ in range(3): f.readline()
# pemapreads = int(f.readline().split()[0])
# pemaprate = f.readline().split()[5].replace("(","")
# L.append(pemaprate)
# f.close()

duprate = round(dupreads/pemapreads,2) # 错误，不应该除以pemapreads，应该除以total

if os.path.isfile(finsertsize):
    f = open(finsertsize)
    for _ in range(7): f.readline()
    meanis = str(int(float(f.readline().split()[5])))
    L.append(meanis)
    f.close()
else:
    L.append("-")

L.append(duprate)

if os.path.isfile(fdepth):
    f = open(fdepth)
    for line in f:
        a,b = line.rstrip().split()
        L.append(a)
        L.append(b)
    f.close()
else:
    L.append("-")
    L.append("-")

#sample, map rate, PE map rate, meanis, duprate, cov1pct, cov10pct
print(sample)
for i in L:
    print(i)

"""output:
mapping rate
PE mapping rate
mean insertsize
dup rate
avg depth
% genome coverage (euchromatic) ≥ 1x  
% genome coverage (euchromatic) ≥ 10x
% genome coverage (euchromatic) ≥ 20x  
% genome coverage (euchromatic) ≥ 30x 
"""
