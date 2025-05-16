import sys, re

sample, splitLog, lfr = sys.argv[1:]

f = open(splitLog)
for line in f:
    if line.startswith('Reads_pair_num ='):
        rp = line.strip().split()[-1]
    elif line.startswith('Reads_pair_num(after split) '):
        match = re.search(r'(\d+)\((\d+\.\d+)%\)', line.strip())
        rp2 = match.group(1)
        bcsplitrate = round(float(match.group(2)),1)
        break
f.close()

D = {}
f = open(lfr)
for line in f:
    k, v = line.strip().split()
    D[k] = v
f.close()
fb, lfrnum, avglen = D['(valid_lfr/valid_barcode)'], D['good-lfr'], D['Average-lfr-length']

print(f"""sample\t{sample}
#read pair\t{rp}
#read pair after split barcode\t{rp2}
barcode split rate\t{bcsplitrate}
fragment/barcode\t{fb}
#fragment (>10k)\t{lfrnum}
fragment avg len\t{avglen}""")