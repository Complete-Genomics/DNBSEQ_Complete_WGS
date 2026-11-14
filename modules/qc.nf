process qc {
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
        
    input:
    val(lib)
    tuple val(id), path(reads) //

    output:
    tuple val(id), path("*qc_1.fq.gz"), path("*qc_2.fq.gz"), emit: reads
    tuple val(id), path("*bssq"), emit: bssq
    tag "$id"
    publishDir "${params.outdir}/$id/fq/", mode: 'link'
    
    script:
    def r1 = "${reads[0]}"
    def r2 = "${reads[1]}"
    """
    mv ${r1} _input_r1.fq.gz
    mv ${r2} _input_r2.fq.gz
    ${params.BIN}SOAPnuke filter \\
      -l 10 -q 0.1 -n 0.01 -T ${task.cpus} \\
      -f CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGA -r TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGG \\
      -1 _input_r1.fq.gz  \\
      -2 _input_r2.fq.gz \\
      -o . \\
      -C ${id}.${lib}.qc_1.fq.gz \\
      -D ${id}.${lib}.qc_2.fq.gz

    mv Basic_Statistics_of_Sequencing_Quality.txt ${id}.${lib}.bssq
    """
    stub:
    "touch ${id}.${lib}.qc_1.fq.gz ${id}.${lib}.qc_2.fq.gz ${id}.${lib}.bssq"
}

process qc_stlfr_stats {
    
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    
    input:
    tuple val(id), path(r1), path(r2) //

    output:
    tuple val(id), env(basecount), emit: basecnt
    tuple val(id), env(rlen), emit: rlen

    tag "$id"
    script:
    """
    ${params.BIN}SOAPnuke filter \\
      -l 0 -q 1  -T ${task.cpus} \\
      -f CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGA -r TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGG \\
      -1 ${r1}  \\
      -2 ${r2} \\
      -o . \\
      -C 1qc.fq.gz \\
      -D 2qc.fq.gz

    rm *qc.fq.gz
    basecount=`head -5 Basic_Statistics_of_Sequencing_Quality.txt | tail -1 | awk '{print \$7}'`
    rlen=`head -2 Basic_Statistics_of_Sequencing_Quality.txt | tail -1 | awk '{print \$3}' | cut -d "." -f 1`
    """
    stub:
    "basecount=1;rlen=100"
}
process readNum {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    
    input:
    tuple val(id), path(bssq) //

    output:
    tuple val(id), env(splitfilenum)

    tag "$id"
    script:
    """
    readnum=`head -3 $bssq | tail -1 | awk '{print \$5}'`
    splitfilenum=`echo "scale=0; \$readnum/${params.splitFqNum}"|bc`
    """
}

