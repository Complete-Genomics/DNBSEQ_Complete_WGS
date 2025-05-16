process frag1 {
    
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), val(chr), path(bam)

    output:
    tuple val(id), path("${id}.stlfr.${chr}.frag1.txt")

    tag "$id, $chr"
    // publishDir "${params.outdir}/$id/"

    script:
    def bam = bam.first()
    """
    # demo.stlfr.chr10.bam -> demo.stlfr.frag1.txt
    sh ${params.SCRIPT}/stat/eachstat_fragment_1.sh $bam . $id 300000 5000 ${params.BIN}samtools 
    mv *.frag1.txt ${id}.stlfr.${chr}.frag1.txt
    """
}

process frag2 {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(frag1s)

    output:
    path "${id}.frag_cov.pdf"
    path "${id}.fraglen_distribution_min5000.pdf"
    path "${id}.frag_per_barcode.pdf"
    tuple val(id), path("*txt"), emit: bcinfo

    tag "$id"

    publishDir "${params.outdir}/report/$id/", mode: 'link'
    // publishDir "${params.outdir}/report/", pattern: "${id}.frag_per_barcode.pdf", saveAs: {"${id}.32.align.fragbc.dist.pdf"}
    // publishDir "${params.outdir}/report/", pattern: "${id}.fraglen_distribution_min5000.pdf", saveAs: {"${id}.33.align.fraglen.dist.pdf"}
    // publishDir "${params.outdir}/report/", pattern: "${id}.frag_cov.pdf", saveAs: {"${id}.34.align.fragcov.dist.pdf"}


    """
    # 
    sh ${params.SCRIPT}/stat/eachstat_fragment_2.sh . $id 5000 ${params.BIN}Rscript
    
    """

}
