import sys, subprocess

prefix, fai = sys.argv[1:]
vcf = prefix + ".phased.vcf.gz"
cmd = f"""
bcftools view -f PASS -v snps -i 'GT="0/1" || GT="1|0" || GT="0|1"' {vcf} |grep -v \# |wc -l
bcftools view -f PASS -v snps -i 'GT="1|0" || GT="0|1"' {vcf} |grep -v \# |wc -l
bcftools view -f PASS -v indels -i 'GT="0/1" || GT="1|0" || GT="0|1"' {vcf} |grep -v \# |wc -l
bcftools view -f PASS -v indels -i 'GT="1|0" || GT="0|1"' {vcf} |grep -v \# |wc -l
tail -3 {prefix}.hapcut_stat.txt |head -1 |awk '{{print $2}}'
grep BLOCK {prefix}.hapblock |wc -l
"""

process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

stdout, stderr = process.communicate()

if process.returncode != 0:
    print(f"Error running command:\n{stderr.decode().strip()}")
else:
    # Get the number of heterozygous SNPs from the command output
    counts = stdout.decode().strip().replace("\n","\t")
    print(counts)







sys.exit()

hetsnp = int(hetsnp)
gsize = int(float(gsize))

f = open(fstat)
for line in f:
    line = line.rstrip()
    if "switch rate" in line:
        switch = line.split()[-1]
    if "phased count" in line:
        phasedcnt = line.split()[-1]
    if line.startswith("N50"):
        n50 = line.split()[-1]
f.close()

hetsnpphased = str(int(phasedcnt)/hetsnp)

spans = 0
f = open(fhb)
for line in f:
    if not line.startswith("BLOCK"):
        continue
    span = line.split()[8]
    spans += int(span)
f.close()

genomephased = spans/gsize

print("\t".join(["Sample","Combined PCR-free + stLFR"]))
print("\t".join(["Phased block N50 (Mb)",n50]))
print("\t".join(["Number of phased blocks",fb]))
print("\t".join(["Switch error rate",switch]))
print("\t".join(["% genome phased",str(genomephased)]))
print("\t".join(["% heterozygous SNPs phased",hetsnpphased]))
