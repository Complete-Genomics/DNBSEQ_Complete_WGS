import sys,os,re,gzip, csv
from funcs import *

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
		cmrg_pct, cmrg_het, cmrg_hom = cmrg_genes(vcf)
		str = f"""
			Sample\t{id}
			Percent of genome coverage >20X (merged bam)\t{merge_genome_cov20}
			Total SNPs called\t{snps:,}
			Total heterozygous SNPs called\t{hetsnps:,}
			Total heterozygous SNPs phased\t{hetsnpsphased:,}
			Total Indels (<50 bp) called\t{indels:,}
			Total heterozygous Indels (<50 bp) called\t{hetindels:,}
			Total phased heterozygous indels\t{hetindelsphased:,}
			Ti/Tv\t{tt}
			Het/hom\t{hh}
			Total cWGS fragments\t{lfrcnt:,}
			Average cWGS length (kb)\t{lfravglen:,}
			Phased contig N50\t{n50:,}
			Total bases in phase block\t{phaseblockbases:,}
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
		_, _, merge_genome_cov20 = fdepth(depthreport)
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
			Phased contig N50\t{n50}
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
			Average fragment length (kb)\t{lfravglen}
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