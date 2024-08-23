process cnv {
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(bam), path(vcf), path(hapblocks)

    output:
    tuple val(id), path("${id}.CNV.result.xls")

    publishDir "${params.outdir}/$id/", mode: 'link'

    script:
    bam = bam.first()
    def vcf = vcf.first()
    def para = (params.ref == "hs37d5") ? "-chr N" : ""
    def a = (params.ref == "hg38" ) ? "-ref GRCH38" : ""
    def pname = (params.var_tool.contains("dv")) ? "${id}.bwa.dv.XXX.hapblock" : "${id}.bwa.gatk.XXX.hapblock"
    """
    ${params.BIN}LFR-cnv -ncpu ${task.cpus} $a  -bam $bam -vcf $vcf -phase . -pname $pname -sp 0.001 -out tmp -lcnv 1000 $para

    mv tmp/ALL.200.format.cnv.1000.highconfidence ${id}.CNV.result.xls
    """

}
process smoove {
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path("out/${id}-smoove.genotyped.vcf.gz")

    tag "$id"
    publishDir "${params.outdir}/$id/", mode: 'link'

    script:
    bam = bam.first()
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    def exclude = "${params.DB}/smoove/exclude_${params.ref}.bed"

    """
    ${params.BIN}smoove call \
      -x --name $id \
      --exclude $exclude \
      --fasta $ref \
      -p ${task.cpus} \
      --genotype $bam \
      --outdir out

    
    """

}

process svpre {
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), val(chr), path(eachbam), path(hapblock)

    output:
    // tuple val(id), val(chr), path("${chr}.vcf"), path("${chr}.region"), path("${chr}.barcode.phase")
    tuple val(id), path("${chr}.*")

    tag "$id, $chr"
    publishDir "${params.outdir}/$id/svsplit/"

    script:
    eachbam = eachbam.first()
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    """
    hapblock=`ls *.hapblock`
    /usr/local/app/stLFRsv/tools/gen_phase/format_phase \
        \$hapblock                 \
        ${chr}.region                       \
        ${chr}.vcf

    /usr/local/app/stLFRsv/tools/gen_phase/get_barcode_from_phase   \
        $ref                                            \
        $eachbam                                        \
        ${chr}.vcf                                      \
        ${chr}.barcode.phase
    """

}
process lfrsv {
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), val(readlen), path(bam), path(svsplits)

    output:
    tuple val(id), path("outdir/*NoRegionFilter"), emit: vcf
    tuple val(id), path("outdir/*bam.final"), emit: bamfinal

    tag "$id"
    publishDir "${params.outdir}/$id/"

    script:
    bam = bam.first()
    // def svsplitdir = "${params.outdir}/$id/svsplit/"
    def svsplitdir = "."
    // def bl = file("${params.DB}/${params.ref}/sv/*bllist")[0]
    // def cl = file("${params.DB}/${params.ref}/sv/*conlist")[0]
    def bl = "${params.DB}/stLFRsv/${params.ref}.bllist"
    def cl = "${params.DB}/stLFRsv/${params.ref}.conlist"
    """
    ${params.BIN}LFR-sv \
      -bam $bam \
      -rlen $readlen \
      -out outdir \
      -ncpu ${task.cpus}  \
      -bl $bl -cl $cl -human Y \
      -phase $svsplitdir

    """

}

process mergesv {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(smoovevcf), path(lfrsvvcf)

    output:
    tuple val(id), path("${id}.SV.vcf")

    tag "$id"
    publishDir "${params.outdir}/$id/"

    script:
    """
    perl ${params.SCRIPT}/stLFRsv_merge_result.pl $lfrsvvcf $smoovevcf 3000 0.7 10000 $id 1 

    head -500 merge.dels.vcf | grep ^# > ${id}.SV.vcf
    cat merge.*.vcf | grep -v ^# | sort -k 1,1 -k 2,2n >> ${id}.SV.vcf
    """

}
