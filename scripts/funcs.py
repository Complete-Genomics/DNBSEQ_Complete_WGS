import sys,os,re,gzip, csv
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
    import pandas as pd
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
    import pandas as pd
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
    import pandas as pd
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
    import pandas as pd
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
