import sys, os, re, datetime


fin, fout = sys.argv[1:]   # id    fqdir:fqdir   1,3-5,22,7-13,20,15-17:1,3-5,22,7-13,20,15-17   hg19

g2 = open("log", 'wt')
def parsebcs(bcs):
	bcs = bcs.split(",")
	for bc in bcs:
		if "~" in bc or "-" in bc:
			a, b = re.split("~|-", bc)
			a, b = int(a), int(b)
			for i in range(a,b + 1):
				bcs.append(str(i))
			bcs.remove(bc)
	return bcs
 
def findfqs(fqdirs, bcs):
	fqdirs = fqdirs.split(",")
	bcs = parsebcs(bcs)
	fq1s, fq2s, fqs = [], [], []
	for bc in bcs:
		for fqdir in fqdirs:
			for f in os.listdir(fqdir):
				if not f.endswith("fq.gz"):
					continue
				suf = re.split(r"_L0[0-9]_", f)[-1]
				bc1 = re.split(r"_[1,2].fq.gz", suf)[0]
				if f.endswith("_1.fq.gz") and (bc1 == bc or bc1 == "0"):
					fq1s.append(os.path.join(fqdir, f))
				elif f.endswith("_2.fq.gz") and (bc1 == bc or bc1 == "0"):
					fq2s.append(os.path.join(fqdir, f))		
	print("findfqs: ", fq1s, fq2s)
	return fq1s, fq2s, fqs
def merge(outName, fqDir, bcs):
	merge = False
	if "," in bcs or "~" in bcs or "-" in bcs:
		merge = True
	fqdirs = [i for i in re.split(':|;',fqDir) if i]
	bcs = [i for i in re.split(':|;',bcs) if i]
	g2.write('\t'.join([str(fqdirs), str(bcs)]) + "\n")
	## merge
	if len(fqdirs) != 1 or len(bcs) != 1:
		merge = True
	if merge:	
		g2.write('merge true' + "\n")
		if len(fqdirs) == 1:	# merge fqs in one fqDir
			bcs = ",".join(bcs)
			fq1s, fq2s, fqs = findfqs(fqdirs[0], bcs)
		else:	# multiple fqDirs
			if len(fqdirs) != len(bcs):
				raise Exception("ERROR: wrong fqdir and barcode info format!")
			fq1s, fq2s, fqs = [], [], []
			for idx in range(len(fqdirs)):
				fqdir = fqdirs[idx]
				bcs1 = bcs[idx]
				a, b, c = findfqs(fqdir, bcs1)
				fq1s += a
				fq2s += b
				fqs += c
		print(fq1s, fq2s, fqs)
		if fq2s and fqs:
			raise Exception("ERROR: mixed SE and PE in specified fqdirs!")
		se = True if fqs else False
		fq1s = fqs if se else fq1s
		if not fq1s:
			raise Exception("ERROR: cant find fastq file with specified fqdir and barcode info")

		s = ".fq.gz" if se else "_1.fq.gz"
		fq1s = " ".join(fq1s)
		outfq1 = outName + s
		command = ' '.join(["cat", fq1s, ">", outfq1])
		t =str(datetime.datetime.now())
		g2.write('\t'.join([t, "merging fq1...", command]) + "\n")
		os.system(command)
		t =str(datetime.datetime.now())
		g2.write('\t'.join([t, "done"])+ "\n")
		if not se:
			s = "_2.fq.gz"
			fq2s = " ".join(fq2s)
			outfq2 = outName + s
			command = ' '.join(["cat", fq2s, ">", outfq2])
			t =str(datetime.datetime.now())
			g2.write('\t'.join([t, "merging fq2...", command])+ "\n")
			os.system(command)
			t =str(datetime.datetime.now())
			g2.write('\t'.join([t, "done"])+ "\n")
		else:
			outfq2 = "-"
		return [outfq1, outfq2]
			
	else:			#no merge; link fq to outdir
		g2.write('merge false' + "\n")
		bc = bcs[0]
		fqdir = fqdirs[0]
		fq1, fq2, fq = findfqs(fqdir, bc)
		if fq:
			return [fq[0], '-']
		elif fq1 and fq2:
			fq1, fq2 = fq1[0], fq2[0]
			command = ["ln -s", fq1, outName + "_1.fq.gz"]
			os.system(' '.join(command))

			command2 = ["ln -s", fq2, outName + "_2.fq.gz"]
			os.system(' '.join(command2))
		else:
			raise Exception("ERROR: cant find fastq file with barcode " + bc)
		return [outName + "_1.fq.gz", outName + "_2.fq.gz"]

	

f = open(fin)
g = open(fout, 'wt')
g.write("sample\tstlfr1\tstlfr2\tpcrfree1\tpcrfree2\n")
for line in f:
	line = line.rstrip().split()
	if "#" in line[0]: continue
	id, fqDir, bcs, ref = line
	g2.write("\t".join([id, fqDir, bcs, ref]) + "\n")

    # items: [[outName_bc1, fq1, fq2],[outName_bc2, fq1, fq2], ...]
    # len(items) equals to len(barcodes) if not merge, otherwise 2
	fq1, fq2 = merge(id, fqDir, bcs)
	fq1, fq2 = os.path.abspath(fq1), os.path.abspath(fq2)
	g.write("\t".join([id, fq1, fq2, fq1, fq2]) + "\n")
f.close()
g.close()
g2.close()
