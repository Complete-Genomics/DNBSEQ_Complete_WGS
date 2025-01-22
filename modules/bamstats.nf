process bamdepth {   
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(lib)
    tuple val(id), path(bam)

    output:
    tuple val(id), env(bamcov)

    tag "$id, $lib"

    script:
    def bam = bam.first()
    "bamcov=`${params.BIN}samtools depth -@ ${task.cpus} $bam | awk '{sum += \$3}END{print sum/${params.ref_len}}'`"
    stub:
    "bamcov=30"
}

process barcodeStat {    
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

    tag "$id"

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path("${id}.barcodeStat.txt") 
 
    script:
    def bam = bam.first()
    """
    ${params.BIN}samtools view $bam|cut -f1,3|sed 's/.*#//g'|awk -v FS="\\t" -v OFS="\\t" '{print \$2,\$1}'| sort --parallel=${task.cpus} | uniq -c |gzip > ${id}.barcodeStat.txt
    """
}
process sameChrBCratio {    
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

    tag "$id"

    input:
    tuple val(id), path(stat)

    output:
    tuple val(id), path("${id}.barcodeStat.txt") 

    // publishDir "${params.outdir}/$id/align/", mode: 'link'
 
    script:
    def bam = bam.first()
    """
    ${params.BIN}Rscript 
    """
}

process insertsize {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(lib)
    tuple val(id), path(bam)

    output:
    path "*pdf"
    tuple val(id), path("*.metrics.txt"), emit: insertsize

    tag "$id, $lib"
    publishDir "${params.outdir}/report/$id/", mode: 'link'

    // publishDir (
    //   path: "${params.outdir}/report/", 
    //   pattern: "*pdf",
    //   saveAs: { fn ->
    //     if (fn.contains("stlfr") && !fn.contains("merge")) { "${id}.43.stlfr.insertsize.dist.pdf" }
    //     else if (fn.contains("pcrfree")) { "${id}.42.pf.insertsize.dist.pdf" }
    //   }
    // )

    script:
    bam = bam.first()
    """    
    ${params.BIN}java -Xms${task.memory.giga}g -Xmx${task.memory.giga}g -jar ${params.SCRIPT}/picard/picard.jar CollectInsertSizeMetrics \\
      I=$bam O=${id}.Insertsize.metrics.txt H=${id}.Insertsize.pdf TMP_DIR=.

    """
    stub:
    "touch ${id}.Insertsize.metrics.txt ${id}.Insertsize.pdf"
}
process gcbias {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(lib)
    tuple val(id), path(bam)

    output:
    path "*pdf"
    tuple val(id), path("*.metrics.txt"), emit: gcbias

    tag "$id, $lib"
    publishDir "${params.outdir}/report/$id/"

    // publishDir (
    //   path: "${params.outdir}/report/", 
    //   pattern: "*pdf",
    //   saveAs: { fn ->
    //     if (fn.contains("stlfr") && !fn.contains("merge")) { "${id}.43.stlfr.insertsize.dist.pdf" }
    //     else if (fn.contains("pcrfree")) { "${id}.42.pf.insertsize.dist.pdf" }
    //   }
    // )

    script:
    bam = bam.first()
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    prefix = "${bam.getBaseName()}"
    """    
    ${params.BIN}java -Xms${task.memory.giga}g -Xmx${task.memory.giga}g -jar ${params.SCRIPT}/picard/picard.jar CollectGcBiasMetrics \\
        R=$ref \\
        I=$bam \\
        O=${id}.GCbias.metrics.txt \\
        CHART=${id}.GCbias.pdf \\
        S=${id}.GCbias.summary_metrics.txt VALIDATION_STRINGENCY=SILENT
    """
    stub:
    "touch ${id}.GCbias.summary_metrics.txt"
}

process samtools_flagstat {
    cpus params.cpu2
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(lib)
    tuple val(id), path(bam)

    output:
    tuple val(id), path("*samtools_flagstat")

    tag "$id, $lib"
    // publishDir "${params.outdir}/$id/align/stats/"

    script:
    def bam = bam.first()
    prefix = "${bam.getBaseName()}"
    "${params.BIN}samtools flagstat -@ ${task.cpus} $bam > ${prefix}.samtools_flagstat"
    stub:
    "touch ${bam.getBaseName()}.samtools_flagstat"
}

process samtools_stats {   
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(lib)
    tuple val(id), path(bam)

    output:
    tuple val(id), path("*.samtools_stats")

    tag "$id, $lib"
    // publishDir "${params.outdir}/$id/align/stats/"

    script:
    def bam = bam.first()
    prefix = "${bam.getBaseName()}"
    """
    ${params.BIN}samtools stats $bam > ${prefix}.samtools_stats
    """
    stub:
    "touch ${bam.getBaseName()}.samtools_stats"
}

process bam2depth {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(lib)
    tuple val(id), path(bam)

    output:
    tuple val(id), path("${id}.${lib}.bam.depth")

    tag "$id, $lib"
    // publishDir "${params.outdir}/$id/align/stats/"

    script:
    def bam = bam.first()
    def non = "${params.DB}/reference/${params.ref}/${params.ref}.nonN.region"
    """
    ${params.BIN}samtools depth -q 0 -Q 0 $bam -b $non > ${id}.${lib}.bam.depth
    """
}

process parseBam2depth { 
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(lib)
    tuple val(id), path(depth)

    output:
    tuple val(id), path("${id}.${lib}.bam.depth.report"), emit: depthreport

    tag "$id, $lib"
    // publishDir "${params.outdir}/$id/align/stats/"

    script:
    def bam = bam.first()
    def non = "${params.DB}/reference/${params.ref}/${params.ref}.nonN.region"
    prefix = "${bam.getBaseName()}"

    """
    python3 ${params.SCRIPT}/calcdepth.py $depth ${params.ref_len} > ${id}.${lib}.bam.depth.report
    """
}
process samtools_depth {
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(lib)
    tuple val(id), path(bam)

    output:
    tuple val(id), path("${id}.${lib}.bamdepth.report")

    tag "$id, $lib"
    // publishDir "${params.outdir}/$id/align/stats/"

    script:
    def bam = bam.first()
    def non = "${params.DB}/${params.ref}/reference/${params.ref}.nonN.region"
    """
    ${params.BIN}samtools depth -@ ${task.cpus} -a -b $non $bam | \\
    awk '\$3 >= 10 { sum10 += 1 } \$3 >= 1 { sum1 += 1 } END { print sum1/${params.ref_len}, sum10/${params.ref_len} }' > ${id}.${lib}.bamdepth.report
    #python3 ${params.SCRIPT}/calcdepth.py ${id}.${lib}.bam.depth ${params.ref_len}> ${id}.${lib}.bamdepth.report
    """
    stub:
    "touch ${id}.${lib}.bamdepth.report"
}

process samtools_depth0 {
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(lib)
    tuple val(id), path(bam)

    output:
    tuple val(id), path("${id}.${lib}.bam.depth")

    tag "$id, $lib"
    // publishDir "${params.outdir}/$id/align/stats/"

    script:
    def bam = bam.first()
    def non = "${params.DB}/reference/${params.ref}/${params.ref}.nonN.region"
    """
    ${params.BIN}samtools depth -@ ${task.cpus} -a -b $non $bam > ${id}.${lib}.bam.depth
    """
}

process splitDepth {
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(depth) //${id}.${lib}.bam.depth

    output:
    tuple val(id), path("*.part") //id.lib.bam.bamdepth.part.xxx

    tag "$id"

    script:
    """
    split -n 50 $depth ${depth.getBaseName()}.bamdepth.part. --additional-suffix=.part
    """
}

process parseDepth {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(depthpart) //demo2.pf.bam.bamdepth.part.bv.part

    output:
    tuple val(id), path("*covbase")

    tag "$id"
    // publishDir "${params.outdir}/$id/align/stats/"

    script:
    def prefix = depthpart.getBaseName() //demo2.pf.bam.bamdepth.part.bv
    """
    awk '\$3 >= 10 { sum10 += \$3 } \$3 >= 1 { sum1 += \$3 } END { print sum1, sum10 }' $depthpart > ${prefix}.covbase
    #python3 ${params.SCRIPT}/calcdepth0.py $depthpart > ${prefix}.report
    """
}

process genomeDepth {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(reports) //demo2.pf.bam.bamdepth.part.bv.report

    output:
    tuple val(id), path("*.depth.report")

    tag "$id"
    // publishDir "${params.outdir}/$id/align/stats/"

    script:
    """
    cat $reports > tmp
    covbase1=`awk '{sum1+=\$1} END{printf "%.0f", sum1}' tmp`
    covbase10=`awk '{sum1+=\$2} END{printf "%.0f", sum1}' tmp`

    cov1=`echo "scale=2; \$covbase1/${params.ref_len}" | bc`
    cov10=`echo "scale=2; \$covbase10/${params.ref_len}" | bc`

    echo \$cov1 > ${reports.first().getBaseName(4)}.depth.report
    echo \$cov10 >> ${reports.first().getBaseName(4)}.depth.report
    rm tmp
    """
}
process eachstat_cov {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(bam) //demo2.pf.bam.bamdepth.part.bv.report

    output:
    tuple val(id), path("*info_1.xls"), emit: info1

    tag "$id"
    publishDir "${params.outdir}/report/$id/"
    // cache false
    script:
    def bam = bam.first()
    def nonN = "${params.DB}/${params.ref}/reference/${params.ref}.nonN.region"
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    """
    perl ${params.SCRIPT}/stat/eachstat_cov.pl $bam ${id}.sorted.bam.info_1.xls samtools
    """
}
process eachstat_depth {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(bam) //demo2.pf.bam.bamdepth.part.bv.report

    output:
    path("*.pdf")
    path("*.xls")
    tuple val(id), path("*info_2.xls"), emit: info2

    tag "$id"
    publishDir "${params.outdir}/report/$id/"
    // cache false
    script:
    def bam = bam.first()
    def nonN = "${params.DB}/${params.ref}/reference/${params.ref}.nonN.region"
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    """
    # output: id.Sequencing.depth.pdf, id.Sequencing.depth.accumulation.pdf, id.sorted.bam.info_2.xls
    perl ${params.SCRIPT}/stat/eachstat_depth.pl $bam . -hg $ref -bd samtools -n $nonN -R . -s $id -t ${task.cpus} > ${id}.sorted.bam.info_2.xls
    """
}
process eachstat_aligncat {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(flgstat), path(stats), path(info1), path(info2), path(insertsize) //demo2.pf.bam.bamdepth.part.bv.report

    output:
    tuple val(id), path("${id}.aligntable.xls")

    tag "$id"
    publishDir "${params.outdir}/report/$id/", mode: 'link'

    script:
    """
    # input: id.sorted.bam.flagstat, id.sorted.bam.stats, id.sorted.bam.info_1.xls, id.sorted.bam.info_2.xls, id.Insertsize.metrics.txt
    perl ${params.SCRIPT}/stat/eachstat_aligncat.pl $flgstat $stats $info1 $info2 $insertsize ${id}.aligntable.xls $id
    """
}
process align_cat {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(lib)
    tuple val(id), path(flagstat), path(stats), path(depth), path(insertsize)

    // flagstat:
    // 1737882 + 0 in total (QC-passed reads + QC-failed reads)
    // 204 + 0 secondary
    // 4996 + 0 supplementary
    // 479 + 0 duplicates
    // 1737158 + 0 mapped (99.96% : N/A)
    // 1732682 + 0 paired in sequencing
    // 866806 + 0 read1
    // 865876 + 0 read2
    // 1700317 + 0 properly paired (98.13% : N/A)
    // 1726960 + 0 with itself and mate mapped
    // 4998 + 0 singletons (0.29% : N/A)
    // 15892 + 0 with mate mapped to a different chr
    // 10125 + 0 with mate mapped to a different chr (mapQ>=5)

    output:
    tuple val(id), path("*.align_cat")

    tag "$id, $lib"
    // publishDir "${params.outdir}/$id/align/stats/"

    script:
    def non = "${params.DB}/reference/${params.ref}/${params.ref}.nonN.region"
    prefix = "${flagstat.getBaseName()}"
    def s
    if (prefix.contains("stlfr") && prefix.contains("merge")) { s = "stLFR"}
    else if (prefix.contains("stlfr")) {s = "stLFR"}
    else {s = "Combined"}
    """
    python3 ${params.SCRIPT}/aligncat.py $flagstat $stats $depth $insertsize $s > ${prefix}.align_cat

    """
    stub:
    "touch ${flagstat.getBaseName()}.align_cat"
}

process align_cat3 {
    
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(lib)
    tuple val(id), path(flagstat), path(stats)

    output:
    tuple val(id), path("*.align_cat")

    tag "$id, $lib"
    // cache false
    publishDir "${params.outdir}/$id/align/stats/"

    script:
    def non = "${params.DB}/reference/${params.ref}/${params.ref}.nonN.region"
    prefix = "${flagstat.getBaseName()}"
    """
    python3 ${params.SCRIPT}/aligncat.py $flagstat $stats fdepth finsertsize Combined > ${prefix}.align_cat
    """
}

process align_catAll {
    
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(align_cat), path(align_cat2), path(align_cat3)

    output:
    path("${id}.align_catall")

    tag "$id"
    // publishDir "${params.outdir}/report/", saveAs: {"${id}.41.align.stats.xls"}

    script:
    """
    echo -e "Library\\nMapping rate\\nPE mapping rate\\nMean insertsize\\nDuplicate rate\\nAverage depth\\n% genome coverage (euchromatic) ≥ 1x\\n% genome coverage (euchromatic) ≥ 10x\\n% genome coverage (euchromatic) ≥ 20x\\n% genome coverage (euchromatic) ≥ 30x\\n" > tmp

    paste tmp $align_cat $align_cat2 $align_cat3 > ${id}.align_catall 
    """
}
