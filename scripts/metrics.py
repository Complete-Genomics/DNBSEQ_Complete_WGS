import sys
from report import bam, cmrg

sample, aligncatstlfr, aligncatpf, fgenecov, stlfrbamdepth, pfbamdepth = sys.argv[1:]

stlfrpemaprate, stlfr_genome_cov10 = bam(aligncatstlfr)
pfpemaprate, pf_genome_cov10 = bam(aligncatpf)
cmrg_cov_pf, cmrg_cov_merge, cmrg_depth_pf, cmrg_depth_merge = cmrg(fgenecov)

str = f"""
    Sample\t{sample}
    stLFR bam avg depth\t{stlfrbamdepth}
    PCR-free bam avg depth\t{pfbamdepth}
    stLFR PE map rate\t{stlfrpemaprate}
    PCR-free PE map rate\t{pfpemaprate}
    stLFR %genome cov > 10x\t{stlfr_genome_cov10}
    PCR-free %genome cov > 10x\t{pf_genome_cov10}
    CMRG avg coverage (PCR-free)\t{cmrg_cov_pf}
    CMRG avg coverage (merged)\t{cmrg_cov_merge}
    CMRG avg depth (PCR-free)\t{cmrg_depth_pf}
    CMRG avg depth (merged)\t{cmrg_depth_merge}
    """

str = '\n'.join(line.strip() for line in str.strip().split('\n'))
print(str)