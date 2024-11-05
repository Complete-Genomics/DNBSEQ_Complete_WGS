import gzip, os, json, subprocess,sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

fin1 = "simCov01.merge.dels.vcf.gz"
fin2 = "simCov05.merge.dels.vcf.gz"

bench = "/hwfssz8/MGI_BIOINFO/USER/biancai/db/HG002_SVs_Tier1_v0.6.vcf.gz"
hg002_smoove = "/hwfssz8/MGI_BIOINFO/USER/biancai/sample4/hg002/smoove_hs37d5/sample-smoove.genotyped.vcf.gz"
hg002_smoove = "/hwfssz8/MGI_BIOINFO/USER/biancai/sample4/hg002/smoove_hs37d5/sample-smoove.genotyped.fixed.vcf.gz"
hg002_svimasm = "/hwfssz9/MGI_LATVIA/BIT/eamozheiko/out/svmasm.bancai/sorted_input.vcf.gz"
hg002_svimasm = "hg002_svimasm.vcf.gz" #from $data/reannotateChr.sh
cmrgsv = "HG002_GRCh37_CMRG_SV_v1.00.vcf.gz"
ref = "/hwfssz8/MGI_BIOINFO/USER/biancai/stLFR_workflow/stLFR_reSeq_v1.4/db/reference/hs37d5/hs37d5.fa"

truvari1="/hwfssz8/MGI_BIOINFO/USER/biancai/Truvari-1.3.4/truvari/truvari"
truvari3 = "/hwfssz8/MGI_BIOINFO/USER/biancai/pwd/envs/sv/bin/truvari"
truvari4 = "/hwfssz8/MGI_BIOINFO/USER/biancai/pwd/envs/truvari/bin/truvari"
bcftools = "/hwfssz8/MGI_BIOINFO/USER/biancai/pwd/bin/bcftools"
tabix = "/hwfssz8/MGI_BIOINFO/USER/biancai/pwd/bin/tabix"

truvari = truvari4

size_ranges = [
    float('-inf'), -30000, -10000, -7500, -5000, -2500, -1000, -750, -500, -250, -50,
    50, 250, 500, 750, 1000, 2500, 5000, 7500, 10000, 30000, float('inf')
]
size_labels = [
    "<-30k", "-30k~-10k", "-10k~-7.5k", "-7.5k~-5k", "-5k~-2.5k", "-2.5k~-1k",
    "-1k~-750", "-750~-500", "-500~-250", "-250~-50", "-50~50", "50~250", "250~500",
    "500~750", "750~1k", "1k~2.5k", "2.5k~5k", "5k~7.5k", "7.5k~10k", "10k~30k", ">30k"
]

def categorize_sv_length(length):
    for i in range(1, len(size_ranges)):
        if size_ranges[i-1] < length <= size_ranges[i]:
            return size_labels[i-1]
    return None
def svcount(vcf, passflg = 1):
    sv_counts = {label: 0 for label in size_labels}
    f = gzip.open(vcf, 'rt')
    for line in f:
        if line.startswith("#"): continue
        line = line.strip().split()
        if passflg:
            if line[6] != "PASS": continue
        info = line[7]
        svtype, svlen = None, None
        for field in info.split(";"):
            if field.startswith("SVTYPE"):
                svtype = field.split("=")[1]
            if field.startswith("SVLEN"):
                svlen = int(field.split("=")[1])
        if svtype and svlen:
            cat = categorize_sv_length(svlen)
            if cat:
                sv_counts[cat] += 1
    f.close()
    sv_counts["-50~50"] = 0
    counts = [sv_counts[label] for label in size_labels]
    return counts
def TRcount(vcf, passflg = 1):
    sv_counts = {label: 0 for label in size_labels}
    f = gzip.open(vcf, 'rt')
    for line in f:
        if line.startswith("#"): continue
        line = line.strip().split()
        if passflg:
            if line[6] != "PASS": continue
        info = line[7]
        svtype, svlen = None, None
        for field in info.split(";"):
            if field.startswith("SVTYPE"):
                svtype = field.split("=")[1]
            if field.startswith("SVLEN"):
                svlen = int(field.split("=")[1])
            if field.startswith("TRall"): #TRUE if at least 20% of the REF bases are tandem repeats of any length
                tr = field.split("=")[1]
        cat = categorize_sv_length(svlen)
        if tr == "TRUE":
            sv_counts[cat] += 1
    f.close()
    sv_counts["-50~50"] = 0
    counts = [sv_counts[label] for label in size_labels]
    return counts
    
def splitbenchvcf():
    cmd = [bcftools, "view -i 'SVTYPE=\"DEL\"'", bench, "-Oz -o", "hg002_bench_dels.vcf.gz"]
    os.system(" ".join(cmd))
    os.system(" ".join([tabix, "-f hg002_bench_dels.vcf.gz"]))

def runTruvari(vcf):
    sample = vcf.split(".")[0]
    b = "hg002_bench_dels.vcf.gz"
    for i in range(11, 21):
        smin = size_ranges[i]
        smax = 200000 if size_ranges[i + 1] == float('inf') else size_ranges[i + 1]
        smin, smax = str(smin), str(smax) 
        outdir = "truvari_dels_" + sample + "_" + smin + "_" + smax
        os.system("rm -rf " + outdir)
        # $truvari -b $bench -c $vcf -o $outdir --reference $ref --passonly -s=50 --sizemax=200000
        cmd = [truvari1, "-b", b, "-c", vcf, "-o", outdir, "--reference", ref, "--passonly", "-s=" + smin, "--sizemax=" + smax]
        if truvari == truvari4:
            cmd = [truvari, "bench -b", b, "-c", vcf, "-o", outdir, "--reference", ref, "--passonly", "-s=" + smin, "--sizemax=" + smax, "--dup-to-ins"]
        print(" ".join(cmd))
        code = os.system(" ".join(cmd))
        if code:
            break
        # break
def parseJson(fin):
    sample = fin.split(".")[0]
    evalD = {'precision':[], 'recall':[]} 
    prefix = "truvari_dels_" + sample
    for i in range(11, 21):
        smin = size_ranges[i]
        smax = 200000 if size_ranges[i + 1] == float('inf') else size_ranges[i + 1]
        smin, smax = str(smin), str(smax) 
        f = prefix + "_" + smin + "_" + smax + "/summary.json"
        with open(f, 'r') as file:
            data = json.load(file)
            precision, recall, f1 = data["precision"], data["recall"], data["f1"]
            tmp = []
            for i in [precision, recall, f1]:
                a = round(i, 3) if i else None
                tmp.append(a)
            precision, recall, f1 = tmp
            evalD['precision'].append(precision)
            evalD['recall'].append(recall)
            # print(t, smin, smax, precision, recall, f1)
    
    a = evalD["precision"][:10][::-1] + [None] + evalD["precision"][10:]
    b = evalD["recall"][:10][::-1] + [None] + evalD["recall"][10:]
    return {"precision":a, "recall": b}

def plotSV():
    realsvcounts1 = svcount(fin1, passflg=0)
    realsvcounts2 = svcount(fin2, passflg=0)
    recall_precision_data = parseJson(fin1)
    recall_precision_data2 = parseJson(fin2)
    # print(len(benchsvcounts),len(recall_precision_data["precision"]), len(recall_precision_data["recall"]))
    # print(benchsvcounts, realsvcounts)



    fig, ax1 = plt.subplots(figsize=(15, 6))
    
    ax1.bar(size_labels, realsvcounts1, label='sim cov 0.1', color='green', alpha=0.3)
    ax1.bar(size_labels, realsvcounts2, label='sim cov 0.5', color='red', alpha=0.3)

    ax1.set_xlabel('SV Size Range')
    ax1.set_ylabel('SV Count')
    ax1.tick_params(axis='y', labelcolor='black')
    ax1.set_xticklabels(size_labels, rotation=45, ha='right')

    ax2 = ax1.twinx()
    print(recall_precision_data)
    print(recall_precision_data2)
    ax2.plot(size_labels, recall_precision_data['recall'], color='green', marker='o', label='Recall')
    ax2.plot(size_labels, recall_precision_data['precision'], color='red', marker='x', label='Precision')
    ax2.plot(size_labels, recall_precision_data2['recall'], color='green', marker='o', label='Recall2', alpha=0.3)
    ax2.plot(size_labels, recall_precision_data2['precision'], color='red', marker='x', label='Precision2', alpha=0.3)
    ax2.set_ylabel('Recall/Precision', color='black')
    ax2.tick_params(axis='y', labelcolor='black')

    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='upper left')

    # plt.legend()
    # plt.title('SV size range')
    # plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for better readability
    plt.savefig('svbenchlenrange.pdf', bbox_inches='tight')
    plt.close()
def plotTR():
    benchsvcounts = svcount(bench)
    TRcounts = TRcount(bench)
    plt.figure(figsize=(15, 6))
    plt.bar(size_labels, benchsvcounts, label = 'all benchmark sv', alpha=0.4)
    plt.bar(size_labels, TRcounts, label = 'ref bases TR > 20%', alpha=0.4, color = 'orange')

    plt.xlabel('SV Size Range')
    plt.ylabel('SV Count')
    plt.title('bench SV size range')
    plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for better readability
    plt.legend()
    plt.savefig('benchSVTR.pdf', bbox_inches='tight')
    plt.close()
def plotCmrg():
    f = gzip.open(cmrgsv, 'rt')
    for line in f:
        if line.startswith("#"): continue
        line = line.strip().split()
        chrom, pos, id, ref, alt, q, flt, info, fmt, hg002 = line
        svlen = len(alt) - len(ref)

        cat = categorize_sv_length(svlen)
        if cat:
            sv_counts[cat] += 1
    f.close()
    sv_counts["-50~50"] = 0
    counts = [sv_counts[label] for label in size_labels]
    print(sv_counts, counts)
    plt.figure(figsize=(15, 6))
    plt.bar(size_labels, counts, alpha=0.6)
    plt.xlabel('SV Size Range')
    plt.ylabel('SV Count')
    plt.title('CMRG SV size range')
    plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for better readability
    plt.savefig('svcmrglenrange.pdf', bbox_inches='tight')
    plt.close()
def isec():
    for i in range(11, 21):
        smin = size_ranges[i]
        smax = 200000 if size_ranges[i + 1] == float('inf') else size_ranges[i + 1]
        smin, smax = str(smin), str(smax) 
        asm = "truvari_svimasm_dels_" + smin + "_" + smax + "/tp-comp.vcf.gz"
        lumpy = "truvari_dels_" + smin + "_" + smax + "/tp-comp.vcf.gz"
        cmd = [bcftools, "isec -C", asm, lumpy, "|wc -l"]
        result = subprocess.run(" ".join(cmd), shell=True, text=True, capture_output=True)
        cnt = result.stdout.strip()
        print("exclusively in asm: ", smin, smax, cnt)

        cmd = [bcftools, "isec -C", lumpy, asm, "|wc -l"]
        result = subprocess.run(" ".join(cmd), shell=True, text=True, capture_output=True)
        cnt = result.stdout.strip()
        print("exclusively in lumpy: ", smin, smax, cnt)

# isec()
# splitbenchvcf()
# runTruvari(fin1)
# runTruvari(fin2)
plotSV()
# plotCmrg()
# plotTR()