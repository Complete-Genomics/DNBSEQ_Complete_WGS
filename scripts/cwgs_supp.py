
import subprocess
import os
import pandas as pd
from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--module", type=str, required=False)
parser.add_argument("--input_fasta", type=str, required=False)

args = parser.parse_args()
module = args.module
input_fasta = args.input_fasta


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


if __name__ == "__main__":

    if module =='nonN':
        count_nonN(input_fasta)

    else:
        print('module not found')