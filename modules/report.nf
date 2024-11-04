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

    ${params.BIN}python3 ${params.SCRIPT}/report.py $id $aligner $varcaller ${id}.bcftoolsStats.txt het $splitLog $lfr $aligncatstlfr $aligncatpf $phase $genecov $vcfeval $vcfevalPf $stlfrbamdepth $pfbamdepth > ${id}.${aligner}.${varcaller}.report
    """
    stub:
    "touch ${id}.${aligner}.${varcaller}.report"
}
process report01 {
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

    ${params.BIN}python3 ${params.SCRIPT}/report_frombam.py $id $aligner $varcaller ${id}.bcftoolsStats.txt het $aligncatstlfr $aligncatpf $phase $genecov $vcfeval $vcfevalPf $stlfrbamdepth $pfbamdepth > ${id}.${aligner}.${varcaller}.report
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
    // println(reports)
    """
    paste -d "," ${params.SCRIPT}/report_template.txt $reports > report.csv
    # ${params.BIN}python3 ${params.SCRIPT}/paste.py ${params.DB}/report_template.txt $reports > report.csv
    """
    stub:
    "touch report.csv"
}

process html {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(flg), path(flg2)

    output:
    path "*html"

    tag "$id"
    publishDir "${params.outdir}/report/", mode: 'link'
    
    // cache false
    script:
    """
    ${params.BIN}convert -resize 750x750 ${params.outdir}/report/$id/${id}.frag_cov.pdf ${id}.frag_cov.png
    ${params.BIN}convert -resize 750x750 ${params.outdir}/report/$id/${id}.fraglen_distribution_min5000.pdf ${id}.fraglen_distribution_min5000.png
    ${params.BIN}convert -resize 750x750 ${params.outdir}/report/$id/${id}.frag_per_barcode.pdf ${id}.frag_per_barcode.png
    ${params.BIN}convert -resize 750x750 ${params.outdir}/report/$id/${id}.frag_cov.pdf ${id}.frag_cov.png
    ${params.BIN}convert -resize 750x750 ${params.outdir}/report/$id/${id}.GCbias.pdf ${id}.GCbias.png
    ${params.BIN}convert -resize 750x750 ${params.outdir}/report/$id/${id}.Insertsize.pdf ${id}.Insertsize.png
    ${params.BIN}convert -resize 750x750 ${params.outdir}/report/$id/${id}.Sequencing.depth.accumulation.pdf ${id}.Sequencing.depth.accumulation.png
    ${params.BIN}convert -resize 750x750 ${params.outdir}/report/$id/${id}.Sequencing.depth.pdf ${id}.Sequencing.depth.png
    ${params.BIN}convert -resize 750x750 ${params.outdir}/report/$id/${id}.haplotype.pdf ${id}.haplotype.png
    ${params.BIN}convert -resize 750x750 ${params.SCRIPT}/stat/legend_circos.pdf ${id}.legend_circos.png

    ${params.BIN}python ${params.SCRIPT}/html/generate_DNApipe_report.py 'stLFR-reSeq V2.0.0.0' $id . ${params.outdir}/report/$id/ .
    """
}
