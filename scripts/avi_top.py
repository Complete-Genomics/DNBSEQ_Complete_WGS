"""Top-N AlphaGenome variant impact (AVI) variants from a VCF annotated with INFO/AVI.

usage: bcftools view <sample>.avi.vcf.gz | python3 avi_top.py <top_n> <top.tsv> <summary.tsv>

AVI is a PHRED score: 10 ~ top 10% of all possible SNVs, 20 ~ top 1%, 30 ~ top 0.1%.
If the VCF also has VEP INFO/CSQ, the most severe consequence is added to each row.
"""
import sys, heapq

IMPACT_RANK = {'HIGH': 3, 'MODERATE': 2, 'LOW': 1, 'MODIFIER': 0}
CSQ_COLS = ['SYMBOL', 'Consequence', 'IMPACT', 'CLIN_SIG', 'AF']
THRESHOLDS = [10, 20, 30]


def avi_value(info):
    # Number=A: one value per ALT; keep the largest annotated one
    vals = [float(v) for v in info.get('AVI', '.').split(',') if v not in ('.', '')]
    return max(vals) if vals else None


def worst_csq(info, csq_fields):
    if not csq_fields or 'CSQ' not in info:
        return ['.'] * len(CSQ_COLS)
    best = None
    for entry in info['CSQ'].split(','):
        rec = dict(zip(csq_fields, entry.split('|')))
        if best is None or IMPACT_RANK.get(rec.get('IMPACT'), -1) > IMPACT_RANK.get(best.get('IMPACT'), -1):
            best = rec
    return [best.get(c) or '.' for c in CSQ_COLS]


def main():
    top_n, top_out, summary_out = int(sys.argv[1]), sys.argv[2], sys.argv[3]
    csq_fields = None
    n_var = n_snv = n_avi = 0
    n_over = {t: 0 for t in THRESHOLDS}
    heap = []  # (avi, seq, row) min-heap of size top_n

    for line in sys.stdin:
        if line.startswith('##INFO=<ID=CSQ'):
            csq_fields = line.split('Format: ')[1].split('"')[0].strip().split('|')
            continue
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        n_var += 1
        ref, alts = f[3], f[4].split(',')
        if len(ref) == 1 and all(len(a) == 1 and a != '*' for a in alts):
            n_snv += 1
        info = dict(kv.split('=', 1) if '=' in kv else (kv, '') for kv in f[7].split(';'))
        avi = avi_value(info)
        if avi is None:
            continue
        n_avi += 1
        for t in THRESHOLDS:
            if avi >= t:
                n_over[t] += 1
        gt = f[9].split(':')[0] if len(f) > 9 else '.'
        row = [f[0], f[1], ref, f[4], f'{avi:g}', gt] + worst_csq(info, csq_fields)
        item = (avi, n_var, row)
        if len(heap) < top_n:
            heapq.heappush(heap, item)
        elif avi > heap[0][0]:
            heapq.heapreplace(heap, item)

    with open(top_out, 'w') as out:
        out.write('\t'.join(['Rank', 'CHROM', 'POS', 'REF', 'ALT', 'AVI', 'GT'] + CSQ_COLS) + '\n')
        for i, (_, _, row) in enumerate(sorted(heap, key=lambda x: (-x[0], x[1])), 1):
            out.write('\t'.join([str(i)] + row) + '\n')

    with open(summary_out, 'w') as out:
        out.write('Metric\tValue\n')
        out.write(f'Variants in VCF\t{n_var:,}\n')
        out.write(f'SNVs\t{n_snv:,}\n')
        out.write(f'Variants with AVI\t{n_avi:,}\n')
        for t in THRESHOLDS:
            out.write(f'AVI >= {t}\t{n_over[t]:,}\n')


if __name__ == '__main__':
    main()
