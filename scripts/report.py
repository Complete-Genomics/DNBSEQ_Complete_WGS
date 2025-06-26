import sys,os,re,gzip, csv
import pandas as pd
from collections import defaultdict
from intervaltree import IntervalTree

def varcnt(varstats):
	*cnt, = open(varstats).readline().strip().split()
	return cnt
def vcfstats(vcf):
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

	return tt, hh

def fvcfeval(vcfeval):
	L = []
	f = open(vcfeval)
	f.readline()
	for line in f:
		for i in line.rstrip().split()[1:]:
			L.append(i)
	f.close()
	return L
def fsplitrate(splitLog):
	f = open(splitLog)
	for line in f:
		if line.startswith('Real_Barcode_types'):
			bctype = line.strip().split()[-1]
		if line.startswith('Reads_pair_num(after split)'):
			bcsplitrate = re.split(r'\(|%\)', line.strip().split()[-1])[1]
			bcsplitrate = f"{round(float(bcsplitrate), 2)}%"
			break
	f.close()
	return bctype, bcsplitrate
def flfr(lfr):
	f = open(lfr)
	for line in f:
		a = line.rstrip().split()[1]
		if "good" in line:
			lfrnum = a
		elif "length" in line:
			avglen = a
		elif "readpair" in line:
			avgfragreadcount = a
		else:
			lfrperbc = a
	f.close()
	avglen = int(float(avglen) / 1e3)
	return lfrnum, avglen
def parse_flgstat(flgstat):
	f = open(flgstat)
	for line in f:
		match = re.search(r'(\d+) \+ \d+ properly paired \(([\d.]+%) : \S+\)', line)
		if match:
			pemaprate = match.group(2)
			break
	f.close()
	return pemaprate

def fphase(phase):
	f = open(phase)
	allhetsnp, phasedhetsnp, allhetindel, phasedhetindel, fbn50, fbnum = f.readline().rstrip().split()
	f.close()
	fbn50 = int(float(fbn50) / 1e6)
	return fbn50

def fblock(hapblock):
	bases = 0
	pattern = re.compile(r'SPAN:\s*(\d+)')

	f = open(hapblock)
	for line in f:
		if not line.startswith('BLOCK'): continue
		match = pattern.search(line)
		bases += int(match.group(1))
	f.close()
	return bases

def fdepth(depthreport):
	cov1, cov10, cov20 = open(depthreport).readline().strip().split()
	cov1 = f"{round(float(cov1) * 100, 1)}%"
	cov10 = f"{round(float(cov10) * 100, 1)}%"
	cov20 = f"{round(float(cov20) * 100, 1)}%"
	return cov1, cov10, cov20

def cmrg(histbed, meanbed):
	ratios = []
	f = open(histbed)
	for line in f:
		line = line.rstrip().split()
		ratios.append(float(line[-1]))
	f.close()
	cov = round(sum(ratios) / len(ratios) * 100, 1)

	depths = []
	f = open(meanbed)
	for line in f:
		line = line.strip().split()
		depths.append(float(line[-1]))
	f.close()
	depth = round(sum(depths) / len(depths) * 100, 1)

	return cov, depth

def load_gene_bed(bed_path):
    """解析基因 BED 文件，返回染色体索引的区间树"""

    # 读取 BED (0-based)
    df = pd.read_csv(bed_path, sep='\t', header=None, 
                     names=['chr', 'start', 'end', 'gene'])
    
    # 构建染色体索引的区间树
    trees = {}
    for _, row in df.iterrows():
        chrom = row['chr']
        if chrom not in trees:
            trees[chrom] = IntervalTree()
        trees[chrom][row['start']:row['end']] = row['gene']  # 左闭右开区间
    return trees

def parse_hapblock(hap_path):
    """解析 hapblock 文件，提取相位变异"""
    blocks = []
    current_block = []
    
    with open(hap_path) as f:
        for line in f:
            if line.startswith("BLOCK"):
                if current_block:  # 保存上一个区块
                    blocks.append(pd.DataFrame(current_block))
                    current_block = []
            elif not line.startswith("********") and line.strip():
                parts = line.split('\t')
                # 提取核心字段：染色体、位置、基因型、质量值
                current_block.append({
                    'chr': parts[3],
                    'pos': int(parts[4]),
                    'gt': parts[7],  # 例如 "0/1:7:3:1,2:0.666667:6,0,22"
                })
    
    return pd.concat(blocks) if blocks else pd.DataFrame()

def classify_variant(gt_str):
    """根据 GT 字段识别纯合/杂合变异"""
    gt = gt_str.split(':')[0]  # 提取 GT 部分 (e.g., "0/1")
    
    if gt in {"1/1", "1|1"}: 
        return 'homozygous'
    elif gt in {"0/1", "1/0", "0|1", "1|0"}: 
        return 'heterozygous'
    return None
def map_variants_to_genes(variants_df, gene_trees):
    """将变异映射到重叠的基因"""
    results = []
    
    for _, var in variants_df.iterrows():
        chrom = var['chr']
        pos = var['pos']
        
        if chrom in gene_trees:
            # 查找包含该位置的基因 (资料 12)
            overlapping_genes = gene_trees[chrom].at(pos)
            for gene_interval in overlapping_genes:
                gene = gene_interval.data
                var_type = classify_variant(var['gt'])
                
                if var_type:
                    results.append({
                        'gene': gene,
                        'pos': pos,
                        'type': var_type,
                        'phase': '|' in var['gt']  # 是否相位已知
                    })
    
    return pd.DataFrame(results)

def count_homozygous_genes(mapped_df):
    """统计含纯合编码变异的基因"""
    return mapped_df[mapped_df['type'] == 'homozygous']['gene'].nunique()
def count_dual_hetero_genes(mapped_df):
    """统计每个等位基因均有杂合变异的基因"""
    
    # 按基因分组相位杂合变异
    gene_alleles = defaultdict(lambda: {'allele0': False, 'allele1': False})
    
    for _, row in mapped_df.iterrows():
        if row['type'] == 'heterozygous' and row['phase']:
            gt = row['gt'].split(':')[0]
            allele0, allele1 = gt.split('|')
            
            # 标记变异所属等位基因
            if allele1 == '1': 
                gene_alleles[row['gene']]['allele1'] = True
            if allele0 == '1': 
                gene_alleles[row['gene']]['allele0'] = True
    
    # 统计双等位基因均有变异的基因
    return sum(1 for alleles in gene_alleles.values() 
              if alleles['allele0'] and alleles['allele1'])

##
def parse_hapblock_blocks(hap_path):
    """解析hapblock文件，返回相位区块的物理边界"""
    blocks = []
    current_block = []
    
    with open(hap_path) as f:
        for line in f:
            if line.startswith("BLOCK"):
                if current_block:  # 处理上一个区块
                    chrom = current_block[0]['chr']
                    positions = [var['pos'] for var in current_block]
                    min_pos = min(positions)
                    max_pos = max(positions)
                    blocks.append({
                        'chrom': chrom,
                        'start0': min_pos - 1,  # 1-based转0-based
                        'end0': max_pos,        # 保持右开区间
                        'span': max_pos - min_pos + 1
                    })
                    current_block = []
            elif not line.startswith("********") and line.strip():
                parts = line.split('\t')
                current_block.append({
                    'chr': parts[3],
                    'pos': int(parts[4])
                })
    return blocks
def build_block_tree(blocks):
    """构建染色体索引的相位区块区间树"""
    block_tree = {}
    for block in blocks:
        chrom = block['chrom']
        if chrom not in block_tree:
            block_tree[chrom] = IntervalTree()
        block_tree[chrom][block['start0']:block['end0']] = block
    return block_tree
def calculate_gene_coverage(gene_trees, block_tree):
    """计算被单个相位区块覆盖的基因"""
    total_genes = 0
    covered_genes = 0
    
    for chrom, tree in gene_trees.items():
        for gene_interval in tree:
            gene_name = gene_interval.data
            gene_start = gene_interval.begin
            gene_end = gene_interval.end
            total_genes += 1
            
            # 检查该染色体是否有相位区块
            if chrom in block_tree:
                # 查询覆盖基因的所有区块
                overlapping_blocks = block_tree[chrom].overlap(gene_start, gene_end)
                covering_blocks = []
                
                for block_interval in overlapping_blocks:
                    block = block_interval.data
                    # 检查是否完全包含基因
                    if block['start0'] <= gene_start and block['end0'] >= gene_end:
                        covering_blocks.append(block)
                
                # 仅被一个区块覆盖
                if len(covering_blocks) == 1:
                    covered_genes += 1
    
    return covered_genes, total_genes

def cmrg_genes():
	gene_trees = load_gene_bed('bed')
	hap_variants = parse_hapblock('hapblock')

	mapped_variants = map_variants_to_genes(hap_variants, gene_trees)
	if mapped_variants.empty:
		cmrg_hom = cmrg_het = 0
	else:
		cmrg_hom = count_homozygous_genes(mapped_variants)
		cmrg_het = count_dual_hetero_genes(mapped_variants)

	##
	blocks = parse_hapblock_blocks("hapblock")
	block_tree = build_block_tree(blocks)
	covered_genes, total_genes = calculate_gene_coverage(gene_trees, block_tree)
	cmrg_pct = f"{round((covered_genes / total_genes) * 100, 1)}%"

	return cmrg_pct, cmrg_het, cmrg_hom

def main():
	flg, *files = sys.argv[1:]
		
	if flg == '0': 
		id, vcf, lfr, histbed, meanbed, depthreport, phasereport = files

		snps, indels, hetsnps, hetindels, hetsnpsphased, hetindelsphased = varcnt('varstat') 

		## lfr
		lfrcnt, lfravglen = flfr(lfr) 

		n50 = fphase(phasereport)
		phaseblockbases = fblock('hapblock')
		_, _, merge_genome_cov20 = fdepth(depthreport)
		cmrg_cov_merge, cmrg_depth_merge = cmrg(histbed, meanbed)

		tt, hh = vcfstats(vcf)
		cmrg_pct, cmrg_het, cmrg_hom = cmrg_genes()
		str = f"""
			Sample\t{id}
			Percent of genome coverage >20X (merged bam){merge_genome_cov20}
			Total SNPs called\t{snps}
			Total heterozygous SNPs called\t{hetsnps}
			Total heterozygous SNPs phased\t{hetsnpsphased}
			Total Indels (<50 bp) called\t{indels}
			Total heterozygous Indels (<50 bp) called\t{hetindels}
			Total phased heterozygous indels\t{hetindelsphased}
			Ti/Tv\t{tt}
			Het/hom\t{hh}
			Total cWGS fragments\t{lfrcnt}
			Average cWGS length (kb)\t{lfravglen}
			Phased contig N50 (Mb)\t{n50}
			Total bases in phase block\t{phaseblockbases}
			Average percent coverage of CMRG genes\t{cmrg_cov_merge}
			Average depth of coverage of CMRG genes\t{cmrg_depth_merge}
			Percent of genes covered by single phased contig\t{cmrg_pct}
			Number of genes with a homozygous coding variant\t{cmrg_hom}
			Number of genes with at least one coding heterozygous variant on each allele\t{cmrg_het}
			"""


	elif flg == 'ref':
		id, vcf, lfr, depthreport, phasereport = files

		snps, indels, hetsnps, hetindels, hetsnpsphased, hetindelsphased = varcnt('varstat')

		lfrcnt, lfravglen = flfr(lfr) 

		n50 = fphase(phasereport)
		phaseblockbases = fblock('hapblock')
		merge_genome_cov20 = fdepth(depthreport)
		tt, hh = vcfstats(vcf)

		str = f"""
			Sample\t{id}
			Percent of genome coverage >20X (merged bam){merge_genome_cov20}
			Total SNPs called\t{snps}
			Total heterozygous SNPs called\t{hetsnps}
			Total heterozygous SNPs phased\t{hetsnpsphased}
			Total Indels (<50 bp) called\t{indels}
			Total heterozygous Indels (<50 bp) called\t{hetindels}
			Total phased heterozygous indels\t{hetindelsphased}
			Ti/Tv\t{tt}
			Het/hom\t{hh}
			Total cWGS fragments\t{lfrcnt}
			Average cWGS length (kb)\t{lfravglen}
			Phased contig N50 (Mb)\t{n50}
			Total bases in phase block\t{phaseblockbases}
			"""
	elif flg.startswith('stlfronly'): # FQC report
		id, vcf, lfr, flgstat, phasereport, stlfrbamdepth, stlfrbamdepthreport = files
		stlfrbamdepth = round(float(stlfrbamdepth), 1)

		snps, indels, hetsnps, hetindels, hetsnpsphased, hetindelsphased = varcnt('varstat')
		lfrcnt, lfravglen = flfr(lfr)
		stlfrpemaprate = parse_flgstat(flgstat)
		_, stlfr_genome_cov10, _ = fdepth(stlfrbamdepthreport)
		
		n50 = fphase(phasereport) 

		str = f"""
			Sample\t{id}
			cWGS bam avg depth\t{stlfrbamdepth}
			Total SNPs called\t{snps}
			Total heterozygous SNPs called\t{hetsnps}
			Total Indels (<50 bp) called\t{indels}
			Total heterozygous Indels (<50 bp) called\t{hetindels}
			Total long fragments\t{lfrcnt}
			Average fragment length\t{lfravglen}
			cWGS mapping rate(unfiltered data)\t{stlfrpemaprate}
			Percent of genome covered by >10X long reads\t{stlfr_genome_cov10}
			Total heterozygous SNPs phased\t{hetsnpsphased}
			Total heterozygous Indels phased\t{hetindelsphased}
			Phased contig N50 (Mb)\t{n50} 
			"""
	elif flg == 'frombam':
		id, aligner, varcaller, varstats, het, aligncatstlfr,aligncatpf, phase, fgenecov, vcfeval, vcfevalpf, stlfrbamdepth, pfbamdepth = files

		snps, indels = varcnt(varstats)
		hetsnps, hetindels = hetvarcnt(het)
		stlfrpemaprate, stlfr_genome_cov10 = bam(aligncatstlfr)
		pfpemaprate, pf_genome_cov10 = bam(aligncatpf)
		hetsnps, hetsnpsphased, hetindels, hetindelsphased, n50, phaseblock = fphase(phase)

		tp, fp, fn, prec, reca, f1, indel_tp, indel_fp, indel_fn, indel_prec, indel_reca, indel_f1 = fvcfeval(vcfeval)
		pf_tp, pf_fp, pf_fn, pf_prec, pf_reca, pf_f1, pf_indel_tp, pf_indel_fp, pf_indel_fn, pf_indel_prec, pf_indel_reca, pf_indel_f1 = fvcfeval(vcfevalpf)

		cmrg_cov_pf, cmrg_cov_merge, cmrg_depth_pf, cmrg_depth_merge = cmrg(fgenecov)

		str = f"""
			Sample\t{id}
			stLFR aligner\t{aligner}
			var caller\t{varcaller}
			stLFR bam avg depth\t{stlfrbamdepth}
			PCR-free bam avg depth\t{pfbamdepth}
			snps\t{snps}
			het snps\t{hetsnps}
			indels\t{indels}
			het indels\t{hetindels}
			##snps eval\t#
			TP\t{tp}
			FP\t{fp}
			FN\t{fn}
			precision\t{prec}
			recall\t{reca}
			f1\t{f1}
			##indels eval\t#
			TP\t{indel_tp}
			FP\t{indel_fp}
			FN\t{indel_fn}
			precision\t{indel_prec}
			recall\t{indel_reca}
			f1\t{indel_f1}
			##PCR-free snps eval\t#
			TP\t{pf_tp}
			FP\t{pf_fp}
			FN\t{pf_fn}
			precision\t{pf_prec}
			recall\t{pf_reca}
			f1\t{pf_f1}
			##PCR-free indels eval\t#
			TP\t{pf_indel_tp}
			FP\t{pf_indel_fp}
			FN\t{pf_indel_fn}
			precision\t{pf_indel_prec}
			recall\t{pf_indel_reca}
			f1\t{pf_indel_f1}
			stLFR PE map rate\t{stlfrpemaprate}
			PCR-free PE map rate\t{pfpemaprate}
			stLFR %genome cov > 10x\t{stlfr_genome_cov10}
			PCR-free %genome cov > 10x\t{pf_genome_cov10}
			het snps phased\t{hetsnpsphased}
			het indels phased\t{hetindelsphased}
			phase block count\t{phaseblock}
			phase block N50\t{n50}
			CMRG avg coverage (PCR-free)\t{cmrg_cov_pf}
			CMRG avg coverage (merged)\t{cmrg_cov_merge}
			CMRG avg depth (PCR-free)\t{cmrg_depth_pf}
			CMRG avg depth (merged)\t{cmrg_depth_merge}
			"""
	elif flg == 'frombam_ref': # no vcfeval, 
		id, aligner, varcaller, varstats, het, aligncatstlfr,aligncatpf, phase, stlfrbamdepth, pfbamdepth = files
		snps, indels = varcnt(varstats)
		hetsnps, hetindels = hetvarcnt(het)
		stlfrpemaprate, stlfr_genome_cov10 = bam(aligncatstlfr)
		pfpemaprate, pf_genome_cov10 = bam(aligncatpf)
		hetsnps, hetsnpsphased, hetindels, hetindelsphased, n50, phaseblock = fphase(phase)

		str = f"""
			Sample\t{id}
			stLFR aligner\t{aligner}
			var caller\t{varcaller}
			stLFR bam avg depth\t{stlfrbamdepth}
			PCR-free bam avg depth\t{pfbamdepth}
			snps\t{snps}
			het snps\t{hetsnps}
			indels\t{indels}
			het indels\t{hetindels}
			stLFR PE map rate\t{stlfrpemaprate}
			PCR-free PE map rate\t{pfpemaprate}
			stLFR %genome cov > 10x\t{stlfr_genome_cov10}
			PCR-free %genome cov > 10x\t{pf_genome_cov10}
			het snps phased\t{hetsnpsphased}
			het indels phased\t{hetindelsphased}
			phase block count\t{phaseblock}
			phase block N50\t{n50}
			"""
	elif flg == 'frombam_ref_PFonly':
		id, varstats, het, aligncatpf, pfbamdepth = files
		snps, indels = varcnt(varstats)
		hetsnps, hetindels = hetvarcnt(het)
		pfpemaprate, pf_genome_cov10 = bam(aligncatpf)

		str = f"""
			Sample\t{id}
			PCR-free bam avg depth\t{pfbamdepth}
			snps\t{snps}
			het snps\t{hetsnps}
			indels\t{indels}
			het indels\t{hetindels}
			PCR-free PE map rate\t{pfpemaprate}
			PCR-free %genome cov > 10x\t{pf_genome_cov10}
			"""
	str = '\n'.join(line.strip() for line in str.strip().split('\n'))
	print(str)

if __name__ == "__main__":
    main()