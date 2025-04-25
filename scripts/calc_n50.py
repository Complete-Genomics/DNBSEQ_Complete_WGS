#!/usr/bin/env python3

from __future__ import print_function
# test
# Author : Peter Edge
# Email  : pedge@eng.ucsd.edu

## edited 02/10/2018 to print 'fraction of SNVs phased and to allow for chr7 to match with '7' in different VCF files

# imports
from collections import defaultdict
import argparse,gzip
import sys
import statistics

desc = '''
Calculate statistics on haplotypes assembled using HapCUT2 or similar tools.
Error rates for an assembled haplotype (specified by -v1 and optionally -h1 arguments)
are computed with respect to a "reference" haplotype (specified by -v2 and optionally -h2 arguments).
All files must contain information for one chromosome only!
To compute aggregate statistics across multiple chromosomes, provide files for
each chromosome/contig as an ordered list, using the same chromosome order between flags.

Note: Triallelic variants are supported, but variants with more than 2 alternative alleles
are currently NOT supported. These variants are ignored. Also, variants where the ref and alt
alleles differ between the test haplotype and reference haplotype are skipped.
'''

def parse_args():

    parser = argparse.ArgumentParser(description=desc)
    # paths to samfiles mapping the same ordered set of RNA reads to different genomes
    parser.add_argument('-v1', '--vcf1', nargs='+', type = str, help='A phased, single sample VCF file to compute haplotype statistics on.')
    parser.add_argument('-v2', '--vcf2', nargs='+', type = str, help='A phased, single sample  VCF file to use as the "ground truth" haplotype.')
    parser.add_argument('-h1', '--haplotype_blocks1', nargs='+', type = str, help='Override the haplotype information in "-v1" with the information in this HapCUT2-format haplotype block file. If this option is used, then the VCF specified with -v1 MUST be the same VCF used with HapCUT2 (--vcf) to produce the haplotype block file!')
    parser.add_argument('-h2', '--haplotype_blocks2', nargs='+', type = str, help='Override the haplotype information in "-v2" with the information in this HapCUT2-format haplotype block file. If this option is used, then the VCF specified with -v2 MUST be the same VCF used with HapCUT2 (--vcf) to produce the haplotype block file!')
    parser.add_argument('-i', '--indels', action="store_true", help='Use this flag to consider indel variants. Default: Indels ignored.',default=False)

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

def parse_hapblock_file(hapblock_file,vcf_file,indels=False):

    snp_ix = 0
    vcf_dict = dict()
    CHROM = None
    with gzip.open(vcf_file,'rt') as infile:
        for line in infile:
            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue

            chrom = el[0]

            if CHROM == None:
                CHROM = chrom
            elif chrom != CHROM:
                print("ERROR: VCFs should contain one chromosome per file")
                print("VCF file:             " + vcf_file)
                exit(1)

            consider = True

            a0 = el[3]
            a1 = el[4]
            a2 = None

            if ',' in a1:
                alt_lst = a1.split(',')
                if len(alt_lst) == 2:
                    a1,a2 = alt_lst
                else:
                    consider = False


            genotype = el[9].split(':')[0]

            ## only consider het snps
            if not (len(genotype) == 3 and genotype[0] in ['0','1','2'] and
                    genotype[1] in ['/','|'] and genotype[2] in ['0','1','2']):
                consider = False

            if genotype[0] == genotype[2]:
                consider = False

            if consider and (not indels) and (('0' in genotype and len(a0) != 1) or
                ('1' in genotype and len(a1) != 1) or ('2' in genotype and len(a2) != 1)):
                consider = False

            genomic_pos = int(el[1])-1
            if consider:
                vcf_dict[snp_ix] = genomic_pos   # start pos
            snp_ix += 1

    blocklist = [] # data will be a list of blocks, each with a list tying SNP indexes to haplotypes

    with open(hapblock_file, 'r') as hbf:
        for line in hbf:
            if len(line) < 3: # empty line
                continue
            if 'BLOCK' in line:
                blocklist.append([])
                continue

            el = line.strip().split('\t')
            if len(el) < 3: # not enough elements to have a haplotype
                continue

            snp_ix = int(el[0])-1

            if snp_ix not in vcf_dict:
                continue

            pos = vcf_dict[snp_ix]

            if len(el) >= 5:

                chrom = el[3]
                if chrom != CHROM and chrom.lstrip('chr') != CHROM.lstrip('chr'):
                    print("ERROR: Chromosome in haplotype block file doesn't match chromosome in VCF")
                    print("Haplotype block file: " + hapblock_file)
                    print("VCF file:             " + vcf_file)
                    exit(1)

                pos2 = int(el[4])-1
                assert(pos == pos2) # if present, genomic position in haplotype block file should match that from VCF

            allele1 = el[1]
            allele2 = el[2]
            ref = el[5]
            alt1 = el[6]
            alt2 = None

            if ',' in alt1:
                alt_lst = alt1.split(',')
                if len(alt_lst) == 2:
                    alt1,alt2 = alt_lst
                else:
                    continue

            blocklist[-1].append((snp_ix, pos, allele1, allele2, ref, alt1, alt2))

    return [b for b in blocklist if len(b) > 1]

def parse_vcf_phase(vcf_file, CHROM, indels = False):

    #block = []
    PS_index = None
    blocks = defaultdict(list)

    with gzip.open(vcf_file, 'rt') as vcf:

        for line in vcf:
            if line[0] == '#':
                continue

            el = line.strip().split('\t')
            if len(el) < 10:
                continue
            if len(el) != 10:
                print("VCF file must be single-sample.")
                exit(1)

            # get the index where the PS information is
            for i,f in enumerate(el[8].split(':')):
                if i == 0:
                    assert(f == 'GT')
                if f == 'PS':
                    if PS_index == None:
                        PS_index = i
                    else:
                        assert(PS_index == i)
                    break

    if PS_index == None:
        print("WARNING: PS flag is missing from VCF. Assuming that all phased variants are in the same phase block.")

    with gzip.open(vcf_file, 'rt') as vcf:

        snp_ix = 0

        for line in vcf:
            if line[0] == '#':
                continue

            el = line.strip().split('\t')
            if len(el) < 10:
                continue
            if len(el) != 10:
                print("VCF file must be single-sample.")
                exit(1)

            consider = True

            phase_data = el[9]

            a0 = el[3]
            a1 = el[4]
            a2 = None

            if ',' in a1:
                alt_lst = a1.split(',')
                if len(alt_lst) == 2:
                    a1,a2 = alt_lst
                else:
                    consider = False

            dat = el[9].split(':')
            genotype = dat[0]

            if not (len(genotype) == 3 and genotype[0] in ['0','1','2'] and
                    genotype[1] in ['|'] and genotype[2] in ['0','1','2']):
                consider = False

            if genotype[0] == genotype[2]:
                consider = False

            if consider and (not indels) and (('0' in genotype and len(a0) != 1) or
                ('1' in genotype and len(a1) != 1) or ('2' in genotype and len(a2) != 1)):
                consider = False

            ps = None
            if PS_index == None:
                ps = 1    # put everything in one block
            elif consider and len(dat) > PS_index:
                ps = dat[PS_index]
                if ps == '.':
                    consider = False

            chrom = el[0]

            if chrom != CHROM and chrom.lstrip('chr') != CHROM.lstrip('chr'):
                print("ERROR: Chromosome in reference haplotype VCF doesn't match chromosome in VCF used for phasing")
                print("reference haplotype VCF: " + vcf_file)
                print("{} != {}".format(CHROM, chrom))
                exit(1)

            pos = int(el[1])-1
            if ps != None and consider and phase_data[1] == '|':
                blocks[ps].append((snp_ix, pos, phase_data[0:1], phase_data[2:3], a0, a1, a2))

            snp_ix += 1

    return [v for k,v in sorted(list(blocks.items())) if len(v) > 1]

# given a VCF file, simply count the number of heterozygous SNPs present.
def count_SNPs(vcf_file,indels=False):

    count = 0

    with gzip.open(vcf_file,'rt') as infile:
        for line in infile:
            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue

            a0 = el[3] #ref
            a1 = el[4] #alt
            a2 = None  

            if ',' in a1:
                alt_lst = a1.split(',')
                if len(alt_lst) == 2:
                    a1,a2 = alt_lst
                else:
                    continue

            genotype = el[9][:3]
            #print ("pre: ",line.strip())

            if not (len(genotype) == 3 and genotype[0] in ['0','1','2'] and
                    genotype[1] in ['/','|'] and genotype[2] in ['0','1','2']):
                continue

            if genotype[0] == genotype[2]:
                continue

            if (not indels) and (('0' in genotype and len(a0) != 1) or
                ('1' in genotype and len(a1) != 1) or ('2' in genotype and len(a2) != 1)):
                continue

            count += 1
    return count

# given a VCF file from a single chromosome return the name of that chromosome
# will almost always be from a single chromosome but don't assume that
def get_ref_name(vcf_file):

    with gzip.open(vcf_file,'rt') as infile:
        for line in infile:
            if line[0] == '#':
                continue
            elif len(line.strip().split('\t')) < 5:
                continue
            return line.strip().split('\t')[0]
    # default
    print("ERROR")
    exit(1)

# combine two dicts
def merge_dicts(d1,d2):
    d3 = d2.copy()
    for k, v in d1.items():
        assert(k not in d3)
        d3[k] = v
    return d3

# the "error_result" abstraction and its overloaded addition operator are handy
# for combining results for the same chromosome across blocks (when the "ground truth"
# is a set of blocks rather than trio), and combining results across different chromosomes
# into genome wide stats
class error_result():
    def __init__(self, ref=None,num_snps=None, AN50_spanlst=None,N50_spanlst=None):

        def create_dict(val,d_type):
            new_dict = defaultdict(d_type)
            if ref != None and val != None:
                new_dict[ref] = val #?
            return new_dict

        self.ref          = set() # set of references in this result (e.g. all chromosomes)
        if ref != None:
            self.ref.add(ref)

        # these are things that can be summed for the same reference,
        # e.g. switch counts for separate blocks are additive

        self.AN50_spanlst   = create_dict(AN50_spanlst,   list)
        self.N50_spanlst    = create_dict(N50_spanlst,    list)

        self.num_snps     = create_dict(num_snps,    int)


    # combine two error rate results
    def __add__(self,other):
        new_err = error_result()

        new_err.ref            = self.ref.union(other.ref)
        new_err.AN50_spanlst   = merge_dicts(self.AN50_spanlst,   other.AN50_spanlst)
        new_err.N50_spanlst    = merge_dicts(self.N50_spanlst,    other.N50_spanlst)
        new_err.num_snps       = merge_dicts(self.num_snps,       other.num_snps)
        return new_err

    def get_num_snps(self):
        return sum(self.num_snps.values())
    
    def get_diff_allele(self):
        return sum(self.diff_allele.values())

    def get_phased_count(self):
        return sum(self.phased_count.values())

    # error rate accessor functions
    def get_switch_rate(self):
        switch_count = self.get_switch_count()
        poss_sw = self.get_poss_sw()
        if poss_sw > 0:
            return float(switch_count)/poss_sw
        else:
            return 0

    def get_mismatch_rate(self):
        mismatch_count = self.get_mismatch_count()
        poss_mm = self.get_poss_mm()

        if poss_mm > 0:
            return float(mismatch_count)/poss_mm
        else:
            return 0


    def get_switch_mismatch_rate(self):
        poss_mm = self.get_poss_mm()

        if poss_mm > 0:
            return float(self.get_switch_count() + self.get_mismatch_count())/poss_mm
        else:
            return 0

    def get_flat_error_rate(self):
        flat_count = self.get_flat_count()
        poss_flat = self.get_poss_flat()
        if poss_flat > 0:
            return float(flat_count)/poss_flat
        else:
            return 0
    def get_blk_len(self):
        N50_spanlst = sum(self.N50_spanlst.values(),[])
        blk_len = sum(N50_spanlst)
        return self.get_switch_count()/blk_len * 10 ** 6

    def get_AN50(self):
        AN50 = 0
        AN50_spanlst = sum(self.AN50_spanlst.values(),[])
        AN50_spanlst.sort(reverse=True)
        phased_sum = 0
        for span, phased in AN50_spanlst:
            phased_sum += phased
            if phased_sum > self.get_num_snps()/2.0:
                AN50 = span
                break
        return AN50

    def get_N50_phased_portion(self):
        N50  = 0
        N50_spanlst = sum(self.N50_spanlst.values(),[])
        N50_spanlst.sort(reverse=True)

        L = sum(N50_spanlst)

        total = 0
        for span in N50_spanlst:
            total += span
            if total > L/2.0:
                N50 = span
                break
        return N50

    def get_median_block_length(self):
        spanlst = sum(self.N50_spanlst.values(),[])
        return statistics.median(spanlst)

#    def get_max_blk_snp_percent(self):
#        snps_in_max_blks = sum(self.maxblk_snps.values())
#        sum_all_snps     = self.get_num_snps()
#
#        if sum_all_snps > 0:
#            return float(snps_in_max_blks) / sum_all_snps
#        else:
#            return 0

    def __str__(self):

        s = ('''
num_snps:           {}
AN50:               {}
N50:                {}
             
            '''.format(self.get_num_snps(),self.get_AN50(), self.get_N50_phased_portion()))

        return s


# compute error rates by using phase data in a VCF as ground truth
# requires VCF to have trio phase information
def hapblock_vcf_error_rate_multiple(assembly_files, vcf_files, indels):

    err = error_result()
    for assembly_file, vcf_file in zip(assembly_files, vcf_files):
        err += hapblock_vcf_error_rate(assembly_file, vcf_file, indels)

    return err

# compute error rates by using phase data in a VCF as ground truth
# requires VCF to have phase information
def hapblock_vcf_error_rate(assembly_file, vcf_file, indels):

    # parse and get stuff to compute error rates
    CHROM = get_ref_name(vcf_file)
    a_blocklist = parse_hapblock_file(assembly_file,vcf_file,indels)

    err = error_rate_calc(a_blocklist, vcf_file, indels)
    return err

def error_rate_calc(a_blocklist, vcf_file, indels=False, phase_set=None):

    ref_name    = get_ref_name(vcf_file)
    num_snps = count_SNPs(vcf_file,indels)

    AN50_spanlst   = []
    N50_spanlst    = []

    for blk in a_blocklist:

        first_pos  = -1
        last_pos   = -1
        first_SNP  = -1
        last_SNP   = -1
        blk_phased = 0

        for snp_ix, pos, a1, a2, ref_str, alt1_str, alt2_str in blk:

            #print('{}\t{}\t{}\t{}'.format(snp_ix, pos, a1, a2))
            if not (a1 == '-' or (phase_set != None and snp_ix not in phase_set)):


                blk_phased+=1
                if first_pos == -1:
                    first_pos = pos
                    first_SNP = snp_ix
                last_pos = pos
                last_SNP = snp_ix

        blk_total = last_SNP - first_SNP + 1

        AN50_spanlst.append(((last_pos-first_pos)*(float(blk_phased)/blk_total), blk_phased))
        N50_spanlst.append((last_pos-first_pos))


    total_error = error_result(ref=ref_name,num_snps=num_snps,AN50_spanlst=AN50_spanlst,N50_spanlst=N50_spanlst)

    return total_error

if __name__ == '__main__':

    args = parse_args()

    print(hapblock_vcf_error_rate_multiple(args.haplotype_blocks1, args.vcf1, args.indels))
