#!/usr/bin/env python3
"""Print the five CMRG report fields from files already prepared in the current directory.

Expected in cwd (made by cmrg_from_results.sh, same as report0 makes them):
  <id>.cmrg.hist.bed   bedtools coverage        (fraction of each CMRG gene covered)
  <id>.cmrg.mean.bed   bedtools coverage -mean  (mean depth of each CMRG gene)
  cmrg_exon.vcf        phased VCF x cmrg273_exon.bed (bedtools intersect -wb)
  hapblock, bed        HapCUT2 hapblock, CMRG gene coordinates
Any field that cannot be computed is printed as FAIL.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from funcs import cmrg, cmrg_genes

sample = sys.argv[1]
cov, depth = cmrg(f'{sample}.cmrg.hist.bed', f'{sample}.cmrg.mean.bed')
pct, het, hom = cmrg_genes(None)
print(f"""Sample\t{sample}
Average percent coverage of CMRG genes\t{cov}
Average depth of coverage of CMRG genes\t{depth}
Percent of genes covered by single phased contig\t{pct}
Number of genes with a homozygous coding variant\t{hom}
Number of genes with at least one coding heterozygous variant on each allele\t{het}""")
