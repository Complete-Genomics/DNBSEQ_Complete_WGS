"""
convert hapcut2 phased block to .bed file

usage (standalone):
    python hapcutPS2bed.py --input sample.lariat.dv.hapblock --output sample.hapblock.bed

usage (legacy multi-file mode, kept for compatibility):
    python hapcutPS2bed.py --input_dir phasesplit/ --sample_name SAMPLE [--fai ref.fai]
"""
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--input",       type=str, required=False, help="single merged hapblock file (cwgs pipeline output)")
parser.add_argument("--output",      type=str, required=False, help="output .bed file")
parser.add_argument("--input_dir",   type=str, required=False, help="legacy: directory of per-chr hapblock files")
parser.add_argument("--sample_name", type=str, required=False, help="legacy: sample name")
parser.add_argument("--fai",         type=str, required=False, default="", help="legacy: fai for non-human")
args = parser.parse_args()


def hapcutPS2bed(hapblock_file, out_bed):
    blocklist = []
    with open(hapblock_file, 'r') as hbf:
        for line in hbf:
            if len(line) < 3:
                continue
            if 'BLOCK' in line:
                blocklist.append([])
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue
            chrom = el[3]
            pos   = int(el[4]) - 1
            blocklist[-1].append((chrom, pos))

    with open(out_bed, 'w') as f:
        for blk in blocklist:
            if not blk:
                continue
            chrom      = blk[-1][0]
            first_pos  = blk[0][1]
            last_pos   = blk[-1][1]
            f.write(f'{chrom}\t{first_pos}\t{last_pos}\n')


if __name__ == "__main__":
    if args.input and args.output:
        hapcutPS2bed(args.input, args.output)

    elif args.input_dir and args.sample_name:
        _dir        = args.input_dir
        sample_name = args.sample_name
        fai         = args.fai

        subprocess.call(f'mkdir -p {_dir}/bed', shell=True)

        if fai == '':
            for i in list(range(1, 23)) + ['X']:
                hapblock_file = f'{_dir}/{sample_name}.lariat.dv.chr{i}.hapblock'
                out_bed       = f'{_dir}/bed/{sample_name}.lariat.dv.chr{i}.bed'
                try:
                    hapcutPS2bed(hapblock_file, out_bed)
                except FileNotFoundError:
                    if i == 'X':
                        print('XY sample, no chrX phased')
        else:
            with open(fai) as f:
                for line in f:
                    chrom = line.strip().split()[0]
                    try:
                        hapcutPS2bed(
                            f'{_dir}/{sample_name}.lariat.dv.{chrom}.hapblock',
                            f'{_dir}/bed/{sample_name}.lariat.dv.{chrom}.bed'
                        )
                    except FileNotFoundError:
                        pass

        subprocess.call(f'cat {_dir}/bed/*.bed > ./{sample_name}.lariat.dv.bed', shell=True)

    else:
        parser.error("provide either --input + --output, or --input_dir + --sample_name")