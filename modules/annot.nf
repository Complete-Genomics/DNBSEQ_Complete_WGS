process vep {
    
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    
    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path("bam.bed")

    script:
    def bam = bam.first()
    """
    ${params.BIN}bedtools bamtobed -i $bam > bam.bed
    """
}