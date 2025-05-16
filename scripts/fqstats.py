import sys,re

"""input looks like:
[biancai@cngb-xcompute-0-17 G400_ECR6_stLFR-1]$ less 01.filter/G400_ECR6_stLFR-1/Basic_Statistics_of_Sequencing_Quality.txt |cat
Item    raw reads(fq1)  clean reads(fq1)        raw reads(fq2)  clean reads(fq2)
Read length     100.0   100.0   100.0   100.0
Total number of reads   818923261 (100.00%)     773902692 (100.00%)     818923261 (100.00%)     773902692 (100.00%)
Number of filtered reads        45020569 (5.50%)        -       45020569 (5.50%)        -
Total number of bases   81892326100 (100.00%)   77390269200 (100.00%)   81892326100 (100.00%)   77390269200 (100.00%)
Number of filtered bases        4502056900 (5.50%)      -       4502056900 (5.50%)      -
Number of base A        24702718134 (30.16%)    23375041354 (30.20%)    23909077267 (29.20%)    22744606991 (29.39%)
Number of base C        16693760363 (20.39%)    15687353889 (20.27%)    16865833482 (20.60%)    15905891895 (20.55%)
Number of base G        16254262300 (19.85%)    15370883154 (19.86%)    16683082808 (20.37%)    15619317583 (20.18%)
Number of base T        24240981047 (29.60%)    22956990803 (29.66%)    24434101024 (29.84%)    23120452731 (29.88%)
Number of base N        604256 (0.00%)  0 (0.00%)       231519 (0.00%)  0 (0.00%)
Q20 number      80400624689 (98.18%)    76026074719 (98.24%)    80148316572 (97.87%)    75805933332 (97.95%)
Q30 number      76853234212 (93.85%)    72705798284 (93.95%)    75675454882 (92.41%)    71635513038 (92.56%)

"""
sample, bssq = sys.argv[1:]

f = open(bssq)
for line in f:
    if line.startswith('Total number of reads'):
        _, raw, clean, *_ = line.split('\t')
        rawpenum = raw.split()[0]
        cleanpenum = clean.split()[0]

    if line.startswith('Total number of bases'):
        _, raw, clean, *_ = line.split('\t')
        rawbasenum = int(raw.split()[0]) * 2
        cleanbasenum = int(clean.split()[0]) * 2

    pattern = r"\((.*?)\)"
    if line.startswith('Number of base C'):
        c1,_,c2,_ = re.findall(pattern, line)
        c1=float(c1.replace("%",""))
        c2 = float(c2.replace("%",""))
    if line.startswith('Number of base G'):
        g1,_,g2,_ = re.findall(pattern, line)
        g1=float(g1.replace("%",""))
        g2=float(g2.replace("%",""))

    if line.startswith('Q20'):
        c1,_,c2,_ = re.findall(pattern, line)
        c1=float(c1.replace("%",""))
        c2 = float(c2.replace("%",""))
        q20 = round((c1 + c2)/2,1)
    if line.startswith('Q30'):
        c1,_,c2,_ = re.findall(pattern, line)
        c1=float(c1.replace("%",""))
        c2 = float(c2.replace("%",""))
        q30 = round((c1 + c2)/2,1)

    # print('Raw reads %bases ≥ Q30' + "\t" + str(round(q30,2)) + "%")
    # print('Total depth' + "\t" + str(round(rawbasenum/3000000000,2)))
f.close()

cg = (c1 + c2 + g1 + g2)/2

print(f"""sample\t{sample}
#raw read pair\t{rawpenum}
#raw bases\t{rawbasenum}
#clean read pair\t{cleanpenum}
#clean bases\t{cleanbasenum}
Q20(%)\t{q20}%
Q30(%)\t{q30}%""")