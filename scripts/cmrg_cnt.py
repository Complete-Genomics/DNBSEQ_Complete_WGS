import sys,os,re, pysam
from intervaltree import IntervalTree

def load_bed_regions(bed_path):
    bed_regions = {}
    with open(bed_path) as f:
        for line in f:
            chrom, start, end, name = line.strip().split('\t')[:4]
            start, end = int(start), int(end)
            if chrom not in bed_regions:
                bed_regions[chrom] = IntervalTree()
            bed_regions[chrom].addi(start, end, data=name)
    return bed_regions

def count_cmrg_genes(vcf_path, bed_regions):
    genes_with_homozygous = set()
    genes_with_het = set()

    vcf = pysam.VariantFile(vcf_path)
    for variant in vcf:
        chrom = variant.chrom
        pos = variant.pos
        # 检查是否在基因区域内
        if chrom in bed_regions:
            overlaps = bed_regions[chrom].at(pos)
			if not overlaps: continue
			for sample in variant.samples
            for interval in overlapping_genes:
                gene_name = interval.data
                # 检查是否为纯合变异
                if variant.num_hom_alt > 0:
                    genes_with_homozygous.add(gene_name)
    return len(genes_with_homozygous)