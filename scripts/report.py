import sys,os,re


def varcnt(varstats):
	f = open(varstats)
	for line in f:
		if not line.startswith("SN"): continue
		if "number of SNPs" in line: 
			snp = line.rstrip().split()[-1]
		elif "number of indels" in line:
			indel = line.rstrip().split()[-1]
			break
	f.close()
	return snp, indel
def hetvarcnt(het):
	f = open(het)
	hetsnp, hetindel = f.readline().rstrip().split()
	f.close()
	return hetsnp, hetindel
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
			bcsplitrate = re.split('\(|\%\)', line.strip().split()[-1])[1]
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
	return lfrnum, avglen
def bam(aligncat):
	f = open(aligncat)
	for _ in range(2): f.readline()
	pemaprate = f.readline().rstrip()
	for _ in range(2): f.readline()
	cov10 = f.readline().rstrip().split()[-1]
	cov10 = f"{round(float(cov10) * 100, 1)}%"
	f.close()
	return pemaprate, cov10
def fphase(phase):
	f = open(phase)
	allhetsnp, phasedhetsnp, allhetindel, phasedhetindel, fbn50, fbnum = f.readline().rstrip().split()
	f.close()
	return allhetsnp, phasedhetsnp, allhetindel, phasedhetindel, fbn50, fbnum
def cmrg(fgenecov):
	f = open(fgenecov)
	pfcov, mergecov, pfdepth, mergedepth = f.readline().rstrip().split()
	f.close()
	return pfcov, mergecov, pfdepth, mergedepth

def main():
	flg, *files = sys.argv[1:]
		
	if flg == '0': 
		id, aligner, varcaller, varstats, het, splitLog, lfr, aligncatstlfr,aligncatpf, phase, fgenecov, vcfeval, vcfevalpf, stlfrbamdepth, pfbamdepth = files

		snps, indels = varcnt(varstats)
		hetsnps, hetindels = hetvarcnt(het)

		tp, fp, fn, prec, reca, f1, indel_tp, indel_fp, indel_fn, indel_prec, indel_reca, indel_f1 = fvcfeval(vcfeval)
		pf_tp, pf_fp, pf_fn, pf_prec, pf_reca, pf_f1, pf_indel_tp, pf_indel_fp, pf_indel_fn, pf_indel_prec, pf_indel_reca, pf_indel_f1 = fvcfeval(vcfevalpf)

		## barcode split rate (pe reads num/pe reads num after split)
		bcsplitrate = fsplitrate(splitLog)

		## lfr
		lfrcnt, lfravglen = flfr(lfr)

		stlfrpemaprate, stlfr_genome_cov10 = bam(aligncatstlfr)
		pfpemaprate, pf_genome_cov10 = bam(aligncatpf)

		hetsnps, hetsnpsphased, hetindels, hetindelsphased, n50, phaseblock = fphase(phase)

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
			barcode split rate\t{bcsplitrate}
			LFR count\t{lfrcnt}
			LFR avg len\t{lfravglen}
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

		# cmrg_hom, cmrg_het = cmrg_genes
		# met = f"""
		# 	Sample\t{id}
		# 	%genome coverage >20X (merged bam){merge_genome_cov20}
		# 	Total SNPs called\t{snps}
		# 	Total heterozygous SNPs called\t{hetsnps}
		# 	Total heterozygous SNPs phased\t{hetsnpsphased}
		# 	Total Indels (< 50 bp) called\t{indels}
		# 	Total heterozygous Indels (< 50 bp) called\t{hetindels}
		# 	Total phased heterozygous indels\t{hetindelsphased}
		# 	Ti/Tv\t{tt}
		# 	Het/hom\t{hh}
		# 	Total cWGS fragments\t{lfrcnt}
		# 	Average cWGS length (kb)\t{int(float(lfravglen) / 1e3)}
		# 	Phased contig N50 (Mb)\t{int(float(n50) /1e6)}
		# 	Total bases in phase block\t{phaseblockbases}
		# 	Average percent coverage of CMRG genes\t{cmrg_cov_merge}
		# 	Average depth of coverage of CMRG genes\t{cmrg_depth_merge}
		# 	Percent of genes covered by single phased contig\t{cmrg_pct}
		# 	Number of genes with a homozygous coding variant\t{cmrg_hom}
		# 	Number of genes with at least one coding heterozygous variant on each allele\t{cmrg_het}
		# 	"""

		# g = open(id + '_metrics.xls','w')
		# met = '\n'.join(line.strip() for line in str.strip().split('\n'))
		# g.write(met)
		# g.close()

	elif flg == 'ref':
		id, aligner, varcaller, varstats, het, splitLog, lfr, aligncatstlfr,aligncatpf, phase, stlfrbamdepth, pfbamdepth = files

		snps, indels = varcnt(varstats)
		hetsnps, hetindels = hetvarcnt(het)

		## barcode split rate (pe reads num/pe reads num after split)
		bcsplitrate = fsplitrate(splitLog)

		## lfr
		lfrcnt, lfravglen = flfr(lfr)

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
			barcode split rate\t{bcsplitrate}
			LFR count\t{lfrcnt}
			LFR avg len\t{lfravglen}
			stLFR PE map rate\t{stlfrpemaprate}
			PCR-free PE map rate\t{pfpemaprate}
			stLFR %genome cov > 10x\t{stlfr_genome_cov10}
			PCR-free %genome cov > 10x\t{pf_genome_cov10}
			het snps phased\t{hetsnpsphased}
			het indels phased\t{hetindelsphased}
			phase block count\t{phaseblock}
			phase block N50\t{n50}
			"""
	elif flg.startswith('stlfronly'):
		id, aligner, varcaller, varstats, het, splitLog, lfr, aligncatstlfr,phase, vcfeval, stlfrbamdepth = files

		snps, indels = varcnt(varstats)
		hetsnps, hetindels = hetvarcnt(het)

		bcpersample, bcsplitrate = fsplitrate(splitLog)

		## lfr
		lfrcnt, lfravglen = flfr(lfr)
		stlfrpemaprate, stlfr_genome_cov10 = bam(aligncatstlfr)
		
		tp, fp, fn, prec, reca, f1, indel_tp, indel_fp, indel_fn, indel_prec, indel_reca, indel_f1 = fvcfeval(vcfeval)
		hetsnps, hetsnpsphased, hetindels, hetindelsphased, n50, phaseblock = fphase(phase)

		str = f"""
			Sample\t{id}
			stLFR bam avg depth\t{stlfrbamdepth}
			snps\t{snps}
			het snps\t{hetsnps}
			indels\t{indels}
			het indels\t{hetindels}
			barcode split rate\t{bcsplitrate}
			LFR count\t{lfrcnt}
			LFR avg len\t{lfravglen}
			stLFR PE map rate(unfiltered data)\t{stlfrpemaprate}
			stLFR %genome cov > 10x\t{stlfr_genome_cov10}
			het snps phased\t{hetsnpsphased}
			het indels phased\t{hetindelsphased}
			phase block count\t{phaseblock}
			phase block N50\t{n50}
			barcode detected per sample\t{bcpersample}
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

	str = '\n'.join(line.strip() for line in str.strip().split('\n'))
	print(str)

if __name__ == "__main__":
    main()