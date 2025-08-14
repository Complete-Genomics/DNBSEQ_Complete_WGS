process report0 {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    val(aligner)
    val(varcaller)
    tuple val(id), path(vcf), path(lfr), path(histbed), path(meanbed), path(depthreport), path(phase)

    output:
    path "${id}.*report"

    tag "$id, $aligner, $varcaller"
    publishDir "${params.outdir}/report/$id/"
    // cache false

    script:
    vcf = vcf.first()
    cmrg_exon_bed = "${params.SCRIPT}/cmrg273_exon.bed"
    """
    set +u
    source /usr/local/miniconda3/bin/activate /usr/local/miniconda3/envs/six

    snp=`bcftools view -v snps $vcf |grep -v \\# |wc -l`
    indel=`bcftools view -v indels $vcf |grep -v \\# |wc -l`
    hetsnp=`bcftools view -v snps -i 'GT="0/1" || GT="1|0" || GT="0|1"' $vcf |grep -v \\# |wc -l`
    hetindel=`bcftools view -v indels -i 'GT="0/1" || GT="1|0" || GT="0|1"' $vcf |grep -v \\# |wc -l`
    phasedhetsnp=`bcftools view -v snps -i 'GT="1|0" || GT="0|1"' $vcf |grep -v \\# |wc -l`
    phasedhetindel=`bcftools view -v indels -i 'GT="1|0" || GT="0|1"' $vcf |grep -v \\# |wc -l`

    echo -e "\$snp\\t\$indel\\t\$hetsnp\\t\$hetindel\\t\$phasedhetsnp\\t\$phasedhetindel" > varstat

    ln -s ${params.DB}/hg38/GRCh38_CMRG_benchmark_gene_coordinates.bed bed
    ln -s ${params.outdir}/$id/phase/${id}.lariat.dv.hapblock hapblock

    bedtools intersect -a $vcf -b $cmrg_exon_bed -wb > cmrg_exon.vcf
    # chr1    1046551 .       A       G       47.2    PASS    .       GT:GQ:DP:AD:VAF:MID:PL:PS       1/1:39:20:0,20:1:small_model:47,39,0:.  chr1    1046397 1046735       AGRN


    ${params.BIN}python3 ${params.SCRIPT}/report.py 0 $id $vcf $lfr $histbed $meanbed $depthreport $phase > ${id}.${aligner}.${varcaller}.report
    """
}
process reportref {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    val(aligner)
    val(varcaller)
    tuple val(id), path(vcf), path(lfr), path(depthreport), path(phase)

    output:
    path "${id}.*report"

    tag "$id, $aligner, $varcaller"
    // publishDir "${params.outdir}/$id/"
    // cache false

    script:
    vcf = vcf.first()
    """
    set +u
    source /usr/local/miniconda3/bin/activate /usr/local/miniconda3/envs/six

    snp=`bcftools view -v snps $vcf |grep -v \\# |wc -l`
    indel=`bcftools view -v indels $vcf |grep -v \\# |wc -l`
    hetsnp=`bcftools view -v snps -i 'GT="0/1" || GT="1|0" || GT="0|1"' $vcf |grep -v \\# |wc -l`
    hetindel=`bcftools view -v indels -i 'GT="0/1" || GT="1|0" || GT="0|1"' $vcf |grep -v \\# |wc -l`
    phasedhetsnp=`bcftools view -v snps -i 'GT="1|0" || GT="0|1"' $vcf |grep -v \\# |wc -l`
    phasedhetindel=`bcftools view -v indels -i 'GT="1|0" || GT="0|1"' $vcf |grep -v \\# |wc -l`

    echo -e "\$snp\\t\$indel\\t\$hetsnp\\t\$hetindel\\t\$phasedhetsnp\\t\$phasedhetindel" > varstat

    ln -s ${params.outdir}/$id/phase/${id}.lariat.dv.hapblock hapblock

    ${params.BIN}python3 ${params.SCRIPT}/report.py ref $id $vcf $lfr $depthreport $phase > ${id}.${aligner}.${varcaller}.report
    """
    stub:
    "touch ${id}.${aligner}.${varcaller}.report"
}
process report_stlfronly {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    val(aligner)
    val(varcaller)
    tuple val(id), path(vcf), path(lfr), path(flgstat), path(phase), val(stlfrbamdepth), path(stlfrbamdepthreport)

    output:
    path "${id}.*report"

    tag "$id, $aligner, $varcaller"
    // publishDir "${params.outdir}/$id/"
    // cache false

    script:
    vcf = vcf.first()
    """
    snp=`bcftools view -v snps $vcf |grep -v \\# |wc -l`
    indel=`bcftools view -v indels $vcf |grep -v \\# |wc -l`
    hetsnp=`bcftools view -v snps -i 'GT="0/1" || GT="1|0" || GT="0|1"' $vcf |grep -v \\# |wc -l`
    hetindel=`bcftools view -v indels -i 'GT="0/1" || GT="1|0" || GT="0|1"' $vcf |grep -v \\# |wc -l`
    phasedhetsnp=`bcftools view -v snps -i 'GT="1|0" || GT="0|1"' $vcf |grep -v \\# |wc -l`
    phasedhetindel=`bcftools view -v indels -i 'GT="1|0" || GT="0|1"' $vcf |grep -v \\# |wc -l`

    echo -e "\$snp\\t\$indel\\t\$hetsnp\\t\$hetindel\\t\$phasedhetsnp\\t\$phasedhetindel" > varstat

    ln -s ${params.outdir}/$id/phase/${id}.lariat.dv.hapblock hapblock

    ${params.BIN}python3 ${params.SCRIPT}/report.py stlfronly $id $vcf $lfr $flgstat $phase $stlfrbamdepth $stlfrbamdepthreport > ${id}.${aligner}.${varcaller}.report
    """
    stub:
    "touch ${id}.${aligner}.${varcaller}.report"
}
process report_stlfronly_ref {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    val(aligner)
    val(varcaller)
    tuple val(id), path(vcf), path(splitLog), path(lfr), path(aligncatstlfr), path(phase), val(stlfrbamdepth)

    output:
    path "${id}.*report"

    tag "$id, $aligner, $varcaller"
    // publishDir "${params.outdir}/$id/"
    // cache false

    script:
    vcf = vcf.first()
    """
    ${params.BIN}bcftools stats $vcf > ${id}.bcftoolsStats.txt
    hetsnp=`${params.BIN}bcftools view -v snps -g het $vcf |grep -v \\# |wc -l`
    hetindel=`${params.BIN}bcftools view -v indels -g het $vcf |grep -v \\# |wc -l`
    echo -e "\$hetsnp\\t\$hetindel" > het

    ${params.BIN}python3 ${params.SCRIPT}/report.py stlfronly_ref $id $lfr $aligncatstlfr $phase $stlfrbamdepth > ${id}.${aligner}.${varcaller}.report
    """
    stub:
    "touch ${id}.${aligner}.${varcaller}.report"
}
process report01 { // from bam
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    val(aligner)
    val(varcaller)
    tuple val(id), path(vcf), path(aligncatstlfr), path(aligncatpf), path(phase), path(genecov), path(vcfeval), path(vcfevalPf), val(stlfrbamdepth), val(pfbamdepth)

    output:
    path "${id}.*report"

    tag "$id, $aligner, $varcaller"
    // publishDir "${params.outdir}/$id/"
    // cache false

    script:
    vcf = vcf.first()
    """
    ${params.BIN}bcftools stats $vcf > ${id}.bcftoolsStats.txt
    hetsnp=`${params.BIN}bcftools view -v snps -g het $vcf |grep -v \\# |wc -l`
    hetindel=`${params.BIN}bcftools view -v indels -g het $vcf |grep -v \\# |wc -l`
    echo -e "\$hetsnp\\t\$hetindel" > het

    ${params.BIN}python3 ${params.SCRIPT}/report.py frombam $id $aligner $varcaller ${id}.bcftoolsStats.txt het $aligncatstlfr $aligncatpf $phase $genecov $vcfeval $vcfevalPf $stlfrbamdepth $pfbamdepth > ${id}.${aligner}.${varcaller}.report
    """
    stub:
    "touch ${id}.${aligner}.${varcaller}.report"
}
process report_frombam_ref { // from bam
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    val(aligner)
    val(varcaller)
    tuple val(id), path(vcf), path(aligncatstlfr), path(aligncatpf), path(phase), val(stlfrbamdepth), val(pfbamdepth)

    output:
    path "${id}.*report"

    tag "$id, $aligner, $varcaller"
    // publishDir "${params.outdir}/$id/"
    // cache false

    script:
    vcf = vcf.first()
    """
    ${params.BIN}bcftools stats $vcf > ${id}.bcftoolsStats.txt
    hetsnp=`${params.BIN}bcftools view -v snps -g het $vcf |grep -v \\# |wc -l`
    hetindel=`${params.BIN}bcftools view -v indels -g het $vcf |grep -v \\# |wc -l`
    echo -e "\$hetsnp\\t\$hetindel" > het

    ${params.BIN}python3 ${params.SCRIPT}/report.py frombam_ref $id $aligner $varcaller ${id}.bcftoolsStats.txt het $aligncatstlfr $aligncatpf $phase $stlfrbamdepth $pfbamdepth > ${id}.${aligner}.${varcaller}.report
    """
    stub:
    "touch ${id}.${aligner}.${varcaller}.report"
}
process report_frombam_PFonly { // from bam
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(vcf), path(aligncatpf), val(pfbamdepth)

    output:
    path "${id}.*report"

    tag "$id"
    // publishDir "${params.outdir}/$id/"
    // cache false

    script:
    vcf = vcf.first()
    """
    ${params.BIN}bcftools stats $vcf > ${id}.bcftoolsStats.txt
    hetsnp=`${params.BIN}bcftools view -v snps -g het $vcf |grep -v \\# |wc -l`
    hetindel=`${params.BIN}bcftools view -v indels -g het $vcf |grep -v \\# |wc -l`
    echo -e "\$hetsnp\\t\$hetindel" > het

    ${params.BIN}python3 ${params.SCRIPT}/report.py frombam_ref_PFonly $id ${id}.bcftoolsStats.txt het $aligncatpf $pfbamdepth > ${id}.frombam.PFonly.report
    """
}
process report {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    
    input:
    path(reports)

    output:
    path "report.csv"

    tag "final report"
    publishDir "${params.outdir}", mode: 'link'
    
    script:
    """
    paste $reports | awk -F'\\t' '{
        line = \$1  # 提取公共的第一列
        # 遍历后续列，每隔两列提取一个（即每个文件的第二列）
        for (i= 2; i <= NF; i += 2) {
            line = line "\t" \$i
        }
        print line
    }' > report.csv

    """
}
process FQC {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    
    input:
    path(reports)

    output:
    path "report_fqc.csv"

    publishDir "${params.outdir}", mode: 'link'
    
    script:
    """
    python ${params.SCRIPT}/fqc.py $reports > report_fqc.csv
    """
}
process html {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    
    input:
    val signal

    output:
    path "*.html"
    // path "*.pdf"

    publishDir "${params.outdir}", mode: 'link'
    
    // cache false
    script:
    """
    ${params.BIN}python ${params.SCRIPT}/my_html.py ${params.outdir}/report/
    """
}
