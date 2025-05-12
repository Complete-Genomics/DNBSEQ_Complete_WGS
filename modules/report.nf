process report0 {
    executor = 'local'
	container false

    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

    input:
    val(aligner)
    val(varcaller)
    tuple val(id), path(vcf), path(splitLog), path(lfr), path(aligncatstlfr), path(aligncatpf), path(phase), path(genecov), path(vcfeval), path(vcfevalPf), val(stlfrbamdepth), val(pfbamdepth)

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

    ${params.BIN}python3 ${params.SCRIPT}/report.py 0 $id $aligner $varcaller ${id}.bcftoolsStats.txt het $splitLog $lfr $aligncatstlfr $aligncatpf $phase $genecov $vcfeval $vcfevalPf $stlfrbamdepth $pfbamdepth > ${id}.${aligner}.${varcaller}.report
    """
    stub:
    "touch ${id}.${aligner}.${varcaller}.report"
}
process reportref {
    executor = 'local'
	container false

    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

    input:
    val(aligner)
    val(varcaller)
    tuple val(id), path(vcf), path(splitLog), path(lfr), path(aligncatstlfr), path(aligncatpf), path(phase), val(stlfrbamdepth), val(pfbamdepth)

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

    ${params.BIN}python3 ${params.SCRIPT}/report.py ref $id $aligner $varcaller ${id}.bcftoolsStats.txt het $splitLog $lfr $aligncatstlfr $aligncatpf $phase $stlfrbamdepth $pfbamdepth > ${id}.${aligner}.${varcaller}.report
    """
    stub:
    "touch ${id}.${aligner}.${varcaller}.report"
}
process report_stlfronly {
    executor = 'local'
	container false

    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

    input:
    val(aligner)
    val(varcaller)
    tuple val(id), path(vcf), path(splitLog), path(lfr), path(aligncatstlfr), path(phase), path(vcfeval), val(stlfrbamdepth)

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

    ${params.BIN}python3 ${params.SCRIPT}/report.py stlfronly $id $aligner $varcaller ${id}.bcftoolsStats.txt het $splitLog $lfr $aligncatstlfr $phase $vcfeval $stlfrbamdepth > ${id}.${aligner}.${varcaller}.report
    """
    stub:
    "touch ${id}.${aligner}.${varcaller}.report"
}
process report_stlfronly_ref {
    executor = 'local'
	container false

    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

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

    ${params.BIN}python3 ${params.SCRIPT}/report.py stlfronly_ref $id $aligner $varcaller ${id}.bcftoolsStats.txt het $splitLog $lfr $aligncatstlfr $phase $stlfrbamdepth > ${id}.${aligner}.${varcaller}.report
    """
    stub:
    "touch ${id}.${aligner}.${varcaller}.report"
}
process report01 { // from bam
    executor = 'local'
    container false

    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

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
    executor = 'local'
    container false

    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

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
process report {
    executor = 'local'
	container false

    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
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
            line = line "," \$i
        }
        print line
    }' > report.csv

    """
    stub:
    "touch report.csv"
}

process html {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    path(flg)

    output:
    path "*html"

    tag "$id"
    publishDir "${params.outdir}/report/", mode: 'link'
    
    // cache false
    script:
    """
    ${params.BIN}python ${params.SCRIPT}/html/generate_cwgs_report.py 'CompleteWGS V1' $id . ${params.outdir}/report/$id/ .
    """
}
