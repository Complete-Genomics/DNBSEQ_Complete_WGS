process fq {
	executor = 'local'
	container false
	
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
	input:
	tuple val(id), path(stlfr1), path(stlfr2), path(pf1), path(pf2)
	
	output:
	tuple val(id), path("${id}_stlfr_1.fq.gz"),  path("${id}_stlfr_2.fq.gz"),emit: stlfr
	tuple val(id), path("${id}_pcrfree_1.fq.gz"), path("${id}_pcrfree_2.fq.gz"),emit: pf

	tag "$id"
	// publishDir "${params.outdir}/$id/fq/"	//, mode: 'copy'

	script:
	"""
	mv $stlfr1 ${id}_stlfr_1.fq.gz
	mv $stlfr2 ${id}_stlfr_2.fq.gz
	mv $pf1 ${id}_pcrfree_1.fq.gz
	mv $pf2 ${id}_pcrfree_2.fq.gz
	"""
}
process fq1 {
	executor = 'local'
	container false
	
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
	input:
	tuple val(id), path(stlfr1), path(stlfr2), val(pf1), val(pf2)
	
	output:
	tuple val(id), path("${id}_stlfr_1.fq.gz"),  path("${id}_stlfr_2.fq.gz")

	tag "$id"
	publishDir "."	//, mode: 'copy'

	script:
	"""
	mv $stlfr1 ${id}_stlfr_1.fq.gz
	mv $stlfr2 ${id}_stlfr_2.fq.gz
	"""
}
process lineNum {
    executor = 'local'
    container false

    input:
    tuple val(id), path(bssq) //

    output:
    tuple val(id), env(lineNum)

    tag "$id"
    script:
    """
    readnum=`head -3 $bssq | tail -1 | awk '{print \$5}'`
    lineNum=`echo "scale=0; ((\$readnum + ${params.lariatSplitFqNum} - 1)/${params.lariatSplitFqNum})*4"|bc`
    """
}

process splitfq {
	cpus params.lariatSplitFqNum / 2 
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

	input:
	tuple val(id), path(r1), path(r2)

    tag "$id"
	output:
	tuple val(id), path("${id}*_1.part_*.fq.gz"), emit: fq1s
	tuple val(id), path("${id}*_2.part_*.fq.gz"), emit: fq2s

	"""
	seqkit split2 -j ${task.cpus} -o .fq.gz -O . -p ${params.lariatSplitFqNum} -1 $r1 -2 $r2
	"""
}

process readLen {
    executor = 'local'
    container false

    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(bssq) //

    output:
    tuple val(id), env(readlen)

    tag "$id"
    script:
    """
    readlen=`head -2 $bssq | tail -1 | awk '{print \$3}' | cut -d "." -f 1`
    """
}
// (base) [biancai@cngb-xcompute-0-17 hg002]$ cat CWGS_run/work/8f/471611f78d09adb307f35161515e28/hg002.pf.bssq
// Item    raw reads(fq1)  clean reads(fq1)        raw reads(fq2)  clean reads(fq2)
// Read length     150.0   150.0   150.1   150.1
// Total number of reads   1624317594 (100.00%)    1604525097 (100.00%)    1624317594 (100.00%)    1604525097 (100.00%)
// Number of filtered reads        19792497 (1.22%)        -       19792497 (1.22%)        -
// Total number of bases   243647639100 (100.00%)  240678764550 (100.00%)  243863347710 (100.00%)  240889880700 (100.00%)
// Number of filtered bases        2968874550 (1.22%)      -       2968874550 (1.22%)      -
// Number of base A        72370644905 (29.70%)    71512115470 (29.71%)    71894308962 (29.48%)    71066765593 (29.50%)
// Number of base C        49244137016 (20.21%)    48622090163 (20.20%)    49803315541 (20.42%)    49169910600 (20.41%)
// Number of base G        49444921940 (20.29%)    48820075729 (20.28%)    49203590736 (20.18%)    48558294124 (20.16%)
// Number of base T        72567648694 (29.78%)    71713077770 (29.80%)    72943152696 (29.91%)    72093814683 (29.93%)
// Number of base N        20286545 (0.01%)        11405418 (0.00%)        18979775 (0.01%)        1095700 (0.00%)
// Q20 number      240365542735 (98.65%)   237507732297 (98.68%)   233220826956 (95.64%)   230733154953 (95.78%)
// Q30 number      230573011123 (94.63%)   227912264206 (94.70%)   198793819193 (81.52%)   196903734036 (81.74%)

process basecount {
    executor = 'local'
    container false

    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(bssq) //

    output:
    tuple val(id), env(basecount)

    tag "$id"
    script:
    """
    basecount=`head -5 $bssq | tail -1 | awk '{print \$7}'`
    """
}

process samplePfFq {
    cpus params.CPU0
    memory params.MEM2 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), val(basecount), val(rlen), path(read) //${id}_qc_1.fq.gz

    output:
    tuple val(id), path("*sampled*fq.gz")

    script:
    def cov = "${params.PF_fq_cov}"
    """
    targetreadnum=`echo "scale=0; ${params.ref_len}*$cov/($rlen*2)" | bc -l`
    targetreadbase=`echo "scale=0; " \$targetreadnum*$rlen | bc -l`
    echo ${params.ref_len}, $cov, ${rlen}. $basecount -\\> \$targetreadbase \$targetreadnum*$rlen > log
    if [ $basecount -gt \$targetreadbase ];then
        echo sample >> log
        ${params.BIN}seqtk sample -2 -s111 $read \$targetreadnum  | ${params.BIN}seqtk trimfq -L $rlen - | gzip >  ${read.getBaseName(2)}.sampled.cov${cov}.fq.gz
    else 
        ln -s $read ${read.getBaseName(2)}.sampled.cov${cov}.fq.gz
        echo not sample >> log
    fi
    """
}

process sampleStlfrFq {
    cpus params.CPU0
    memory params.MEM2 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), val(basecount), val(rlen), path(read) //${id}_qc_1.fq.gz

    output:
    tuple val(id), path("*sampled*fq.gz")

    script:
    def cov = "${params.stLFR_fq_cov}"
    """
    targetreadnum=`echo "scale=0; ${params.ref_len}*$cov/($rlen*2)" | bc -l`
    targetreadbase=`echo "scale=0; " \$targetreadnum*$rlen | bc -l`
    echo ${params.ref_len}, $cov, ${rlen}. $basecount -\\> \$targetreadbase \$targetreadnum*$rlen > log
    if [ $basecount -gt \$targetreadbase ];then
        echo sample >> log
        ${params.BIN}seqtk sample -2 -s111 $read \$targetreadnum  | gzip >  ${read.getBaseName(2)}.sampled.cov${cov}.fq.gz
    else 
        ln -s $read ${read.getBaseName(2)}.sampled.cov${cov}.fq.gz
        echo not sample >> log
    fi
    """
}
process fqcheck {
  // tag {fastq}
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(lib)
    tuple val(id), path(r1), path(r2) //hg002_stlfr_1.sampled.cov0.fq.gz

    output:
    tuple val(id), path("${id}.${lib}_1.fqcheck"), path("${id}.${lib}_2.fqcheck") 
    
    tag "$id, $lib"
    // publishDir "${params.outdir}/$id/fq/"

    script:
    """
    ${params.BIN_fqcheck}fqcheck33 -r $r1 -c ${id}.${lib}_1.fqcheck
    ${params.BIN_fqcheck}fqcheck33 -r $r2 -c ${id}.${lib}_2.fqcheck
    """
}
process fqdist {
  
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    publishDir "${params.outdir}/report/$id", mode: 'link'
  // publishDir (
  //   path: "${params.outdir}/report/",
  //   saveAs: { fn ->
  //     if (fn.contains("base.png") && fn.contains("stlfr")) { "${id}.24.qc-stlfr.basecomp.png" }
  //     else if (fn.contains("qual.png") && fn.contains("stlfr")) { "${id}.26.qc-stlfr.basequal.png" }
  //     else if (fn.contains("base.png") && fn.contains("pf")) { "${id}.23.qc-pf.basecomp.png" }
  //     else if (fn.contains("qual.png") && fn.contains("pf")) { "${id}.25.qc-pf.basequal.png" }
  //   }
  // )

  input:
  val(lib)
  tuple val(id), path(fqcheck1), path(fqcheck2)

  output:
  path "*.png"

  tag "$id"
  script:
  //hg002_stlfr_1
  // def prefix = fqcheck1.getBaseName().replaceFirst(/_1/, '')
  """
  perl ${params.SCRIPT}/fqcheck/fqcheck_distribute.pl ${fqcheck1} ${fqcheck2} -o ${id}.${lib}.
  """
}
process eachstat_fastq {
    
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(bssq), path(splitlog), path(bcsumAndReadcount)

    output:
    path("*.xls")

    // publishDir "${params.outdir}/report/", saveAs: {"${id}.31.qc-stlfr.stats.xls"}
    publishDir "${params.outdir}/report/$id/", mode: 'link'

    script:
    """
    # eachstat_fastq.pl <sample> <indir> <output> <genomesize>  <alignstatdir>
    # in: Basic_Statistics_of_Sequencing_Quality.txt, split_stat_read1.log, 02.align/id/stat/min5000frag_barcode_summary.txt, 02.align/id/stat/minfrag5000_mean_len_readcount.txt
    # out: id.fastqtable.xls, id.fragtable.xls

    mv $bssq Basic_Statistics_of_Sequencing_Quality.txt 
    ${params.BIN}perl ${params.SCRIPT}/stat/eachstat_fastq.pl $id ${params.ref_len} 
    """
}
