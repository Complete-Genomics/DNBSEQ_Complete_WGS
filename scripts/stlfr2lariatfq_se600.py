#!/usr/bin/env python3
"""
SE600 barcode split: extract 3x10bp barcodes from read sequence,
encode them in read name (#bc1_bc2_bc3), output gDNA-only FASTQ.

Usage:
    python stlfr2lariatfq_se600.py barcode.list input.fq.gz output.fq.gz umi_start

UMI layout (42 bp): bc1(10) + spacer(6) + bc2(10) + spacer(6) + bc3(10)
  umi_start (1-based) = position where UMI begins in the read
    - UMI at read start: umi_start = 1   (gDNA = seq[42:])
    - UMI at read end:   umi_start = read_len - 41

Barcode matching uses 1-mismatch tolerance (matches split_cwgs.pl).

Performance:
  - dnaio C parser for FASTQ I/O
  - xopen auto-uses pigz/igzip for multi-threaded gzip
  - Falls back to gzip if dnaio not installed
"""

import sys

try:
    import dnaio
    HAS_DNAIO = True
except ImportError:
    HAS_DNAIO = False
    import gzip


NUCS = ('A', 'C', 'G', 'T')
BC_LEN = 10
SPACER = 6
UMI_LEN = 42  # bc1(10) + spacer(6) + bc2(10) + spacer(6) + bc3(10)


def load_barcode_hash(barcode_file):
    """Load barcode list and build seq -> id hash with 1-mismatch tolerance.

    barcode.list format: id,seq (comma-separated, e.g. "1,ACGTAACGTA")
    """
    bc_hash = {}
    with open(barcode_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            if len(parts) < 2:
                # try whitespace fallback
                parts = line.split()
                if len(parts) < 2:
                    continue
                bc_seq, bc_id = parts[0], parts[1]
            else:
                bc_id, bc_seq = parts[0], parts[1]

            for pos in range(len(bc_seq)):
                for nuc in NUCS:
                    variant = bc_seq[:pos] + nuc + bc_seq[pos+1:]
                    bc_hash[variant] = bc_id
    return bc_hash


def process_dnaio(input_fq, output_fq, bc_hash, umi_start_0):
    """Fast path: dnaio + xopen (uses pigz/igzip automatically)."""
    bc1_s = umi_start_0
    bc2_s = umi_start_0 + BC_LEN + SPACER
    bc3_s = umi_start_0 + (BC_LEN + SPACER) * 2
    umi_end = umi_start_0 + UMI_LEN

    n_total = 0
    n_matched = 0

    with dnaio.open(input_fq, mode='r') as reader, \
         dnaio.open(output_fq, mode='w', fileformat='fastq') as writer:
        for record in reader:
            n_total += 1
            seq = record.sequence
            qual = record.qualities

            if len(seq) < umi_end:
                continue

            bc1 = seq[bc1_s:bc1_s + BC_LEN]
            bc2 = seq[bc2_s:bc2_s + BC_LEN]
            bc3 = seq[bc3_s:bc3_s + BC_LEN]

            bc1_id = bc_hash.get(bc1)
            bc2_id = bc_hash.get(bc2)
            bc3_id = bc_hash.get(bc3)

            if bc1_id and bc2_id and bc3_id:
                n_matched += 1
                bc_tag = bc1_id + '_' + bc2_id + '_' + bc3_id
            else:
                bc_tag = '0_0_0'

            gdna_seq = seq[:umi_start_0] + seq[umi_end:]
            gdna_qual = qual[:umi_start_0] + qual[umi_end:]

            old_name = record.name.split()[0]
            if old_name.endswith('/1') or old_name.endswith('/2'):
                old_name = old_name[:-2]
            new_name = old_name + '#' + bc_tag + '/1'

            writer.write(dnaio.SequenceRecord(new_name, gdna_seq, gdna_qual))

    return n_total, n_matched


def process_gzip(input_fq, output_fq, bc_hash, umi_start_0):
    """Fallback path: stdlib gzip (no dnaio installed)."""
    bc1_s = umi_start_0
    bc2_s = umi_start_0 + BC_LEN + SPACER
    bc3_s = umi_start_0 + (BC_LEN + SPACER) * 2
    umi_end = umi_start_0 + UMI_LEN

    n_total = 0
    n_matched = 0

    h_in = gzip.open(input_fq, 'rt')
    wh = gzip.open(output_fq, 'wt')
    line_no = 0
    rec_name = rec_seq = rec_plus = rec_qual = None

    for line in h_in:
        line_no += 1
        line = line.rstrip('\n')
        mod = line_no % 4
        if mod == 1:
            rec_name = line
        elif mod == 2:
            rec_seq = line
        elif mod == 3:
            rec_plus = line
        else:
            rec_qual = line
            n_total += 1
            seq = rec_seq
            qual = rec_qual
            if len(seq) >= umi_end:
                bc1 = seq[bc1_s:bc1_s + BC_LEN]
                bc2 = seq[bc2_s:bc2_s + BC_LEN]
                bc3 = seq[bc3_s:bc3_s + BC_LEN]
                bc1_id = bc_hash.get(bc1)
                bc2_id = bc_hash.get(bc2)
                bc3_id = bc_hash.get(bc3)
                if bc1_id and bc2_id and bc3_id:
                    n_matched += 1
                    bc_tag = bc1_id + '_' + bc2_id + '_' + bc3_id
                else:
                    bc_tag = '0_0_0'
                gdna_seq = seq[:umi_start_0] + seq[umi_end:]
                gdna_qual = qual[:umi_start_0] + qual[umi_end:]
                old_name = rec_name.split()[0].lstrip('@')
                if old_name.endswith('/1') or old_name.endswith('/2'):
                    old_name = old_name[:-2]
                new_name = '@' + old_name + '#' + bc_tag + '/1'
                wh.write(new_name + '\n' + gdna_seq + '\n+\n' + gdna_qual + '\n')

    h_in.close()
    wh.close()
    return n_total, n_matched


def main():
    if len(sys.argv) != 5:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    barcode_file = sys.argv[1]
    input_fq = sys.argv[2]
    output_fq = sys.argv[3]
    umi_start_0 = int(sys.argv[4]) - 1  # convert to 0-based

    if umi_start_0 < 0:
        print(f"ERROR: umi_start must be >= 1 (got {sys.argv[4]})", file=sys.stderr)
        sys.exit(1)

    print(f"Loading barcode hash...", file=sys.stderr)
    bc_hash = load_barcode_hash(barcode_file)
    print(f"  {len(bc_hash)} barcode variants loaded", file=sys.stderr)
    print(f"UMI starts at position {umi_start_0 + 1} (1-based), length {UMI_LEN}",
          file=sys.stderr)

    if HAS_DNAIO:
        print("Using dnaio (fast path)", file=sys.stderr)
        n_total, n_matched = process_dnaio(input_fq, output_fq, bc_hash, umi_start_0)
    else:
        print("dnaio not installed, using stdlib gzip (slow path)", file=sys.stderr)
        n_total, n_matched = process_gzip(input_fq, output_fq, bc_hash, umi_start_0)

    pct = 100.0 * n_matched / n_total if n_total > 0 else 0.0
    print(f"Total reads:   {n_total}", file=sys.stderr)
    print(f"Matched reads: {n_matched} ({pct:.2f}%)", file=sys.stderr)


if __name__ == '__main__':
    main()