process barcode_split {
    cpus params.CPU1
    memory params.MEM2 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
        
    input:
    tuple val(id), path(r1), path(r2)

    output:
    tuple val(id), path("${id}_split_1.fq.gz"), path("${id}_split_2.fq.gz"), emit: reads
    tuple val(id), path("split_stat_read1.log"), emit: log

    tag "$id"
    publishDir "${params.outdir}/$id/fq/", mode: 'link'

    script:
    def bcList = "${params.DB}/barcode/barcode.list"
    """
    fl=`zcat $r1 | head -n 1 || true`
    echo \$fl
    if echo \$fl | grep -q \\#; then 
        echo skip split barcode
        mv $r1 ${id}_split_1.fq.gz
        mv $r2 ${id}_split_2.fq.gz

        ${params.BIN}python ${params.SCRIPT}/splitRate.py ${id}_split_1.fq.gz > split_stat_read1.log
    else
        echo split barcode
        tmp=`zcat $r1 | head -n 2 | tail -n 1 |wc -c ||true`
        rlen=`expr \$tmp - 1`
        echo \$rlen
        mv $r1 read_1.fq.gz
        mv $r2 read_2.fq.gz

        s1=\$((\$rlen*2+1))
        s2=\$((\$rlen*2+17))
        s3=\$((\$rlen*2+33))

        bcpos="-I \$s1 10 1 false -I \$s2 10 1 false -I \$s3 10 1 false"
        ${params.BIN}MGI.Lite.GenFastQ -F read_1.fq.gz read_2.fq.gz \\
                            -B ${bcList}    \\
                            \$bcpos       \\
                            --stLFR                     \\
                            -O split_out  \\
                            --logPath ./ 

        out1=`find split_out -name "*_1.fq.gz"`
        out2=`find split_out -name "*_2.fq.gz"`
        log=`find split_out -name "split_stat_read.log"`

        mv \$out1 ${id}_split_1.fq.gz
        mv \$out2 ${id}_split_2.fq.gz
        mv \$log split_stat_read1.log
    fi
    """
    stub:
    "touch ${id}_split_1.fq.gz ${id}_split_2.fq.gz split_stat_read1.log"
}
process get_rlen {
    cpus params.CPU1
    memory params.MEM2 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
        
    input:
    tuple val(id), path(r1), path(r2)

    output:
    tuple val(id), env(rlen)

    tag "$id"

    script:
    """

    tmp=`zcat $r1 |head -2 |tail -1 | wc -c`
    rlen=\$((\$tmp-1)) 

    s1=\$((\$rlen*2+1))
    s2=\$((\$rlen*2+17))
    s3=\$((\$rlen*2+33))

    bcpos="-I \$s1 10 1 false -I \$s2 10 1 false -I \$s3 10 1 false"
    ${params.BIN}MGI.Lite.GenFastQ -F read_1.fq.gz read_2.fq.gz \\
                        -B ${bcList}    \\
                        \$bcpos       \\
                        --stLFR                     \\
                        -O split_out  \\
                        --logPath ./ 
 
    out1=`find split_out -name "*_1.fq.gz"`
    out2=`find split_out -name "*_2.fq.gz"`
    log=`find split_out -name "split_stat_read.log"`

    mv \$out1 ${id}_split_1.fq.gz
    mv \$out2 ${id}_split_2.fq.gz
    mv \$log split_stat_read1.log
    """
}
process barcode_split_stLFRreseq {
    cpus params.CPU1
    memory params.MEM2 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
        
    input:
    tuple val(id), path(r1), path(r2)

    output:
    tuple val(id), path("${id}_split_1.fq.gz"), path("${id}_split_2.fq.gz"), emit: reads
    tuple val(id), path("split_stat_read1.log"), emit: log
    tuple val(id), env(rlen), emit: rlen

    tag "$id"
    publishDir "${params.outdir}/$id/fq/", mode: 'link'

    script:
    def bcList = "${params.DB}/barcode/barcode.list"
    """
    mv $r1 read_1.fq.gz
    mv $r2 read_2.fq.gz

    tmp=`zcat read_1.fq.gz | head -2 | tail -1 | wc -c || true`
    rlen=\$((\$tmp-1))

    s1=\$((\$rlen*2+1))
    s2=\$((\$rlen*2+17))
    s3=\$((\$rlen*2+33))

    bcpos="-I \$s1 10 1 false -I \$s2 10 1 false -I \$s3 10 1 false"
    ${params.BIN}MGI.Lite.GenFastQ -F read_1.fq.gz read_2.fq.gz \\
                        -B ${bcList}    \\
                        \$bcpos       \\
                        --stLFR                     \\
                        -O split_out  \\
                        --logPath ./ 
 
    out1=`find split_out -name "*_1.fq.gz"`
    out2=`find split_out -name "*_2.fq.gz"`
    log=`find split_out -name "split_stat_read.log"`

    mv \$out1 ${id}_split_1.fq.gz
    mv \$out2 ${id}_split_2.fq.gz
    mv \$log split_stat_read1.log
    """
}
process makeLog {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
        
    input:
    tuple val(id), path(r1), path(r2)

    output:
    tuple val(id), path("split_stat_read1.log")

    tag "$id"
    publishDir "${params.outdir}/$id/fq/", mode: 'link'
    // cache false
    script:
    """
    ${params.BIN}python ${params.SCRIPT}/splitRate.py $r1 > split_stat_read1.log
    """
}