
import subprocess
import os
import gzip
import pandas as pd
from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--module", type=str, required=False)
parser.add_argument("--input_fasta", type=str, required=False)
parser.add_argument("--db_path", type=str, required=False)

args = parser.parse_args()
module = args.module
input_fasta = args.input_fasta
db_path = args.db_path

def count_nonN(input_fasta):
    """
    Reads a .fasta file and outputs non-N regions to a .bed file.

    Parameters:
    input_fasta (str): Path to the input .fasta file.
    output_bed (str): Path to the output .bed file.
    """
    output_bed = f'{input_fasta}.nonN.region'
    with open(output_bed, 'w') as bed_file:
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq = str(record.seq)
            chrom = record.id  # Chromosome or sequence name
            start = None  # Start of the non-N region

            for i, base in enumerate(seq):
                if base.upper() != "N" and start is None:
                    start = i  # Start of a non-N region
                elif base.upper() == "N" and start is not None:
                    # End of a non-N region
                    bed_file.write(f"{chrom}\t{start}\t{i}\n")
                    start = None

            # If the sequence ends with a non-N region
            if start is not None:
                bed_file.write(f"{chrom}\t{start}\t{len(seq)}\n")
    return


def correct_hg38_fasta(db_path):
    """
    cretea GCA_000001405.15_GRCh38_no_alt_analysis_set_corrected.fasta from GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

    """
    hg38_fa = f'{db_path}/_RAWDATA_/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz'
    output_fa = f'{db_path}/hg38/panGenome/GCA_000001405.15_GRCh38_no_alt_analysis_set_corrected.fasta'

    chromosomes_to_remove = {'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr1_KI270711v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrUn_GL000213v1', 'chrUn_KI270338v1', 'chrUn_KI270364v1', 'chrUn_KI270371v1', 'chrUn_KI270372v1', 'chrUn_KI270374v1', 'chrUn_KI270375v1', 'chrUn_KI270424v1', 'chrUn_KI270528v1', 'chrUn_KI270581v1', 'chrUn_KI270587v1'}
    # # Now filter the FASTA file
    print(f"\nFiltering {hg38_fa} to create {output_fa}...")

    open_fn = gzip.open if hg38_fa.endswith('.gz') else open
    with open_fn(hg38_fa, 'rt') as infile, open(output_fa, 'w') as outfile:
        write_current = False
        current_chrom = None

        for line in infile:
            if line.startswith('>'):
                # New sequence header
                current_chrom = line[1:].strip().split()[0]  # Get chromosome name

                if current_chrom in chromosomes_to_remove:
                    write_current = False
                    print(f"  Skipping chromosome: {current_chrom}")
                else:
                    write_current = True
                    outfile.write(line)
            elif write_current:
                # Write sequence lines only for kept chromosomes
                outfile.write(line)
    
    print(f"\nDone! Created {output_fa} with only chromosomes present in both .fai files.")
    count_nonN(output_fa)


if __name__ == "__main__":

    if module =='nonN':
        count_nonN(input_fasta)

    elif module == 'corrected_fasta':
        correct_hg38_fasta(db_path)

    else:
        print('module not found')