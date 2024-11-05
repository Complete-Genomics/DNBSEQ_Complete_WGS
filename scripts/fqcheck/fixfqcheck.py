import os,sys

prefix = sys.argv[1]
fin1 = prefix + "_1.fqcheck"
fin2 = prefix + "_2.fqcheck"

f1 = open(fin1)
for line in f1:
    if line.strip().startswith('A'):
        b1 = line.strip().split()[-1]
        break
f1.close()

f2 = open(fin2)
for line in f2:
    if line.strip().startswith('A'):
        b2 = line.strip().split()[-1]
        break
f2.close()

if b1 != b2:
    b1, b2 = int(b1), int(b2)
    fin = fin1 if b2 > b1 else fin2
    f = open(fin)
    L = range(min(b1,b2) + 1, max(b1,b2) + 1)

    command = ['mv',fin,fin + ".bak"]
    os.system(' '.join(command))

    g = open(fin, 'wt')
    for line in f:
        
        if line.strip().startswith('A'):
            line = line.rstrip() + "\t" + "\t".join([str(i) for i in L]) + "\n"
            g.write(line)
        elif line.strip().startswith('Total') or line.strip().startswith('base'):
            line = line.rstrip() + "\t" + "\t".join(["0" for _ in range(len(L))]) + "\n"
            g.write(line)
        else:
            g.write(line)
            continue

    g.close()
    f.close()