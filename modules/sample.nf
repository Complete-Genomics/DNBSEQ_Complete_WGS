process sample_stlfr {
    
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
        
    input:
    tuple val(id), val(rlen), path(read)

    output:
    tuple val(id), val(rlen), path("*sampled*fq.gz")

    publishDir "${params.outdir}/$id/fq/"

    script:
    def cov = "${params.stLFR_sampling_cov}"
    """
	cov=$cov
    stlfrfqlen=${rlen}
	n=\$((${params.ref_len}*\$cov/(\$stlfrfqlen*2)))
	${params.BIN}seqtk sample -2 -s111 ${read} \$n  | gzip > ${read.getBaseName(2)}.sampled0.cov${cov}.fq.gz
    """
}

process sample_pf {
    
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
        
    input:
    tuple val(id), val(rlen), path(read)

    output:
    tuple val(id), val(rlen), path("*sampled*fq.gz")

    publishDir "${params.outdir}/$id/fq/"

    script:
    def cov = "${params.PF_sampling_cov}"
    """
	cov=$cov
	pffqlen=${rlen}
    n=\$((${params.ref_len}*\$cov/(\$pffqlen*2)))
    ${params.BIN}seqtk sample -2 -s111 $read \$n  | ${params.BIN}seqtk trimfq -L \$pffqlen - | gzip >  ${read.getBaseName(2)}.sampled.cov${cov}.fq.gz
    """
}
