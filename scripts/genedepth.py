import sys

genebed, bed = sys.argv[1:]

genes = []
f = open(genebed)
for line in f:
	chr, st, ed = line.rstrip().split()
	genes.append([chr, int(st), int(ed)])	
f.close()

curr = 0
out = []
d = 0
cnt = 0
curr_region = genes[curr]
f = open(bed)
for line in f:
	_chr,_pos, _depth = line.rstrip().split()
	_pos = int(_pos)
	_depth = int(_depth)
	chr, st, ed = curr_region
	print(chr,st,ed,_pos,_pos>=st,_pos<=ed,_depth)
	if chr == _chr and _pos >= st and _pos <= ed:
		cnt += 1
		d += _depth
	elif chr != _chr or _pos > ed:
		avgd = d/cnt if cnt else 0
		out.append([chr,st,ed,avgd])
		curr += 1
		if curr < len(genes):	
			curr_region = genes[curr]
			cnt = 0
			d = 0			
		else:
			break
f.close()	
for i in out:
	chr, st, ed, avgd = i
	#print("\t".join([chr, str(st), str(ed), str(avgd)]))
