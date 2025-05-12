import sys

sample, bcsplit, lfr = sys.argv[1:]

print(f"""sample\t{sample}
#read pair\t{rp}
#read pair after split barcode\t{rp2}
barcode split rate\t{bcsplitrate}
barcode per fragment\t{bpf}
#fragment (>)\t{lfrnum}
fragment avg len\t{avglen}

""")