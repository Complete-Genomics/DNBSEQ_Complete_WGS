import sys,os,re

#python3 ${params.SCRIPT}/report0.py $id $vcf $splitLog $lfr $aligncatstlfr $aligncatpf $phase $genecov > ${id}.report
"""
[biancai@cngb-xcompute-0-17 G400_ECR6_stLFR-1]$ head 01.filter/G400_ECR6_stLFR-1/split_stat_read1.log
Barcode_types = 1536*1536*1536=3623878656
Barcode_types_with_mismatch = 62976*62976*62976=249761340850176(1)
Real_Barcode_types = 33574842
Reads_pair_num = 818923261
Reads_pair_num(after split) = 753672975(92.032186%)
0       65250286        0_0_0
8282618 1       1000_1000_1064
11119045        1       1000_1000_1096

"""
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
	try:
		for _ in range(4): f.readline()
		bcsplitrate = re.split('\(|\%\)', f.readline().split()[-1])[1]
		bcsplitrate = round(float(bcsplitrate), 2)
	except:
		bcsplitrate = f.readline().strip()
	f.close()
	return bcsplitrate
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

elif flg == 'stlfronly':
	id, aligner, varcaller, varstats, het, splitLog, lfr, aligncatstlfr,phase, vcfeval, stlfrbamdepth = files

	snps, indels = varcnt(varstats)
	hetsnps, hetindels = hetvarcnt(het)

	bcsplitrate = fsplitrate(splitLog)

	## lfr
	lfrcnt, lfravglen = flfr(lfr)
	stlfrpemaprate, stlfr_genome_cov10 = bam(aligncatstlfr)
	
	tp, fp, fn, prec, reca, f1, indel_tp, indel_fp, indel_fn, indel_prec, indel_reca, indel_f1 = fvcfeval(vcfeval)
	hetsnps, hetsnpsphased, hetindels, hetindelsphased, n50, phaseblock = fphase(phase)

	str = f"""
		Sample\t{id}
		stLFR aligner\t{aligner}
		var caller\t{varcaller}
		stLFR bam avg depth\t{stlfrbamdepth}
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
		barcode split rate\t{bcsplitrate}
		LFR count\t{lfrcnt}
		LFR avg len\t{lfravglen}
		stLFR PE map rate\t{stlfrpemaprate}
		stLFR %genome cov > 10x\t{stlfr_genome_cov10}
		het snps phased\t{hetsnpsphased}
		het indels phased\t{hetindelsphased}
		phase block count\t{phaseblock}
		phase block N50\t{n50}
		"""
elif flg == 'stlfronly_ref':
	id, aligner, varcaller, varstats, het, splitLog, lfr, aligncatstlfr,phase, stlfrbamdepth = files
	snps, indels = varcnt(varstats)
	hetsnps, hetindels = hetvarcnt(het)

	bcsplitrate = fsplitrate(splitLog)

	lfrcnt, lfravglen = flfr(lfr)
	stlfrpemaprate, stlfr_genome_cov10 = bam(aligncatstlfr)
	hetsnps, hetsnpsphased, hetindels, hetindelsphased, n50, phaseblock = fphase(phase)

	str = f"""
		Sample\t{id}
		stLFR aligner\t{aligner}
		var caller\t{varcaller}
		stLFR bam avg depth\t{stlfrbamdepth}
		snps\t{snps}
		het snps\t{hetsnps}
		indels\t{indels}
		het indels\t{hetindels}
		barcode split rate\t{bcsplitrate}
		LFR count\t{lfrcnt}
		LFR avg len\t{lfravglen}
		stLFR PE map rate\t{stlfrpemaprate}
		stLFR %genome cov > 10x\t{stlfr_genome_cov10}
		het snps phased\t{hetsnpsphased}
		het indels phased\t{hetindelsphased}
		phase block count\t{phaseblock}
		phase block N50\t{n50}
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