
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

    common_chroms = {'chrUn_KI270580v1', 'chr14_GL000225v1_random', 'chrUn_KI270579v1', 'chrUn_KI270521v1', 'chr16_KI270728v1_random', 'chr1_KI270712v1_random', 'chrUn_KI270742v1', 'chrUn_KI270442v1', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270755v1', 'chr7', 'chr1_KI270706v1_random', 'chrUn_GL000219v1', 'chrUn_GL000214v1', 'chrUn_KI270412v1', 'chrUn_KI270529v1', 'chrUn_KI270394v1', 'chrUn_KI270467v1', 'chrUn_KI270391v1', 'chr2_KI270716v1_random', 'chrUn_KI270468v1', 'chrUn_KI270420v1', 'chrUn_KI270517v1', 'chrUn_KI270393v1', 'chrY', 'chr9_KI270720v1_random', 'chr6', 'chrUn_KI270395v1', 'chr2_KI270715v1_random', 'chrUn_KI270448v1', 'chr12', 'chrUn_KI270311v1', 'chr9_KI270717v1_random', 'chrUn_KI270303v1', 'chr14_GL000009v2_random', 'chrUn_KI270378v1', 'chr17_KI270730v1_random', 'chrUn_KI270390v1', 'chrUn_KI270336v1', 'chr5_GL000208v1_random', 'chrUn_KI270417v1', 'chrUn_KI270522v1', 'chr17', 'chrUn_KI270741v1', 'chrUn_KI270530v1', 'chrUn_KI270510v1', 'chr14', 'chrUn_KI270548v1', 'chrUn_KI270516v1', 'chrUn_KI270757v1', 'chrUn_KI270310v1', 'chr3_GL000221v1_random', 'chrUn_KI270392v1', 'chr19', 'chr13', 'chrUn_KI270388v1', 'chrUn_KI270322v1', 'chrUn_KI270418v1', 'chrUn_KI270511v1', 'chr9_KI270719v1_random', 'chrUn_KI270305v1', 'chrX', 'chr11_KI270721v1_random', 'chrUn_KI270429v1', 'chrUn_KI270539v1', 'chrUn_GL000218v1', 'chr17_GL000205v2_random', 'chr2', 'chr16', 'chrUn_GL000195v1', 'chrUn_KI270749v1', 'chrUn_KI270518v1', 'chr1', 'chr1_KI270710v1_random', 'chrUn_KI270538v1', 'chrUn_GL000226v1', 'chrUn_KI270317v1', 'chrUn_KI270414v1', 'chr22_KI270734v1_random', 'chrUn_KI270753v1', 'chrUn_KI270466v1', 'chrUn_KI270381v1', 'chrY_KI270740v1_random', 'chr8', 'chr9', 'chr22_KI270738v1_random', 'chrUn_GL000224v1', 'chrUn_KI270383v1', 'chr14_KI270725v1_random', 'chrUn_KI270411v1', 'chrUn_KI270333v1', 'chrUn_KI270396v1', 'chrUn_KI270582v1', 'chrUn_KI270512v1', 'chrUn_KI270379v1', 'chr1_KI270709v1_random', 'chrUn_KI270363v1', 'chrUn_KI270304v1', 'chrUn_KI270465v1', 'chrUn_KI270382v1', 'chrUn_KI270593v1', 'chr17_KI270729v1_random', 'chr9_KI270718v1_random', 'chrUn_KI270384v1', 'chrUn_KI270337v1', 'chrUn_KI270752v1', 'chrUn_KI270312v1', 'chr3', 'chr14_KI270722v1_random', 'chrUn_KI270386v1', 'chrEBV', 'chrUn_KI270389v1', 'chrUn_KI270507v1', 'chrUn_KI270419v1', 'chrUn_KI270340v1', 'chr5', 'chrUn_KI270362v1', 'chr1_KI270714v1_random', 'chr20', 'chrUn_KI270746v1', 'chrUn_KI270745v1', 'chrUn_KI270435v1', 'chrUn_KI270583v1', 'chrUn_KI270590v1', 'chrUn_KI270335v1', 'chrUn_KI270329v1', 'chrUn_KI270747v1', 'chrUn_KI270425v1', 'chrUn_KI270748v1', 'chrUn_KI270751v1', 'chrUn_KI270519v1', 'chrUn_KI270387v1', 'chrM', 'chr18', 'chrUn_KI270385v1', 'chrUn_KI270754v1', 'chrUn_KI270373v1', 'chr1_KI270713v1_random', 'chrUn_KI270316v1', 'chrUn_KI270509v1', 'chrUn_KI270744v1', 'chr1_KI270708v1_random', 'chr14_KI270724v1_random', 'chr15', 'chrUn_KI270315v1', 'chr11', 'chr21', 'chrUn_KI270750v1', 'chrUn_KI270334v1', 'chr1_KI270707v1_random', 'chrUn_KI270515v1', 'chrUn_KI270302v1', 'chrUn_GL000220v1', 'chrUn_KI270422v1', 'chrUn_KI270588v1', 'chrUn_GL000216v2', 'chrUn_KI270508v1', 'chr4_GL000008v2_random', 'chr22_KI270733v1_random', 'chr10', 'chrUn_KI270584v1', 'chr14_KI270723v1_random', 'chrUn_KI270376v1', 'chrUn_KI270366v1', 'chrUn_KI270756v1', 'chrUn_KI270320v1', 'chrUn_KI270330v1', 'chr22_KI270732v1_random', 'chrUn_KI270544v1', 'chrUn_KI270589v1', 'chr22', 'chrUn_KI270423v1', 'chr22_KI270731v1_random', 'chrUn_KI270743v1', 'chr4', 'chrUn_KI270438v1', 'chrUn_KI270591v1', 'chr14_GL000194v1_random'}
    chromosomes_to_remove= ['chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr1_KI270711v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrUn_GL000213v1', 'chrUn_KI270338v1', 'chrUn_KI270364v1', 'chrUn_KI270371v1', 'chrUn_KI270372v1', 'chrUn_KI270374v1', 'chrUn_KI270375v1', 'chrUn_KI270424v1', 'chrUn_KI270528v1', 'chrUn_KI270581v1', 'chrUn_KI270587v1']
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
                
                if current_chrom in common_chroms:
                    write_current = True
                    outfile.write(line)
                else:
                    write_current = False
                    print(f"  Skipping chromosome: {current_chrom}")
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