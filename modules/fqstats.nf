process fqstats {
    
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(bssq) //id.pcrfree.bssq

    output:
    path("${id}.fqstats.xls")

    tag "$id"
    publishDir "${params.outdir}/report/$id/", mode: 'link'
     
    script:
    """
    python3 ${params.SCRIPT}/fqstats.py $id $bssq > ${id}.fqstats.xls
    """
    stub:
    "touch ${id}.fqstats.xls"
}

// process fqstats_stlfr {
//     
//     clusterOptions = "-clear -cwd -l vf=${params.MEM1}g,num_proc=${params.CPU0}  -binding linear:${params.CPU0} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

//     input:
//     tuple val(id), path(bssq)

//     output:
//     path("out")

//     publishDir "${params.outdir}/report/"//, saveAs: "${id}.22.qc-stlfr.xls"
 
//     script:
//     """
//     python3 $baseDir/../bin/fqstats.py $bssq > out
//     """
// }
