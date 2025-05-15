import sys, gzip
from report import varcnt, hetvarcnt, fphase

sample, vcf, varstats, het, phase = sys.argv[1:]

snps, indels = varcnt(varstats)
hetsnps, hetindels = hetvarcnt(het)
hetsnps, hetsnpsphased, hetindels, hetindelsphased, n50, phaseblock = fphase(phase)
n50 = round(int(n50) / 1000000, 1)

ti = tv = het = hom = 0

f = gzip.open(vcf, 'rt')
for line in f:
    if line.startswith('#'): continue
    line = line.rstrip().split()
    ref, alt, info = line[3], line[4], line[-1]
    gt = info.split(':')[0]

    if len(ref) != 1 or len(alt) != 1: continue
    if ref not in {"A", "T", "C", "G"} or alt not in {"A", "T", "C", "G"}: continue

    if (ref, alt) in {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}:
        ti += 1
    else:
        tv += 1
    
    a1, a2 = gt.replace('|', '/').split('/')
    if a1 == '.' or a2 == '.': continue
    if a1 != a2:
        het += 1
    elif a1 == '1':
        hom += 1
f.close()


tt = ti / tv if tv > 0 else 0
hh = het / hom if hom > 0 else 0

tt = round(tt, 1)
hh = round(hh, 1)

str = f"""
    Sample\t{sample}
    #SNPs\t{snps}
    #het SNPs\t{hetsnps}
    #indels\t{indels}
    #het indels\t{hetindels}
    #phased SNPs\t{hetsnpsphased}
    #phased indels\t{hetindelsphased}
    Tv/Ti\t{tt}
    Het/Hom\t{hh}
    #phase block\t{phaseblock}
    phased contig N50 (Mb)​\t{n50}
    """

str = '\n'.join(line.strip() for line in str.strip().split('\n'))
print(str)