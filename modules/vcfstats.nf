process vcfstats {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(vcf), path(phase)

    output:
    tuple val(id), path("${id}.vcfstats.xls")

    tag "$id"
    publishDir "${params.outdir}/report/$id/"

    script:
    vcf = vcf.first()
    """
    ${params.BIN}bcftools stats $vcf > ${id}.bcftoolsStats.txt
    hetsnp=`${params.BIN}bcftools view -v snps -g het $vcf |grep -v \\# |wc -l`
    hetindel=`${params.BIN}bcftools view -v indels -g het $vcf |grep -v \\# |wc -l`
    echo -e "\$hetsnp\\t\$hetindel" > het

    python3 ${params.SCRIPT}/vcfstats.py $id $vcf ${id}.bcftoolsStats.txt het $phase > ${id}.vcfstats.xls
    """
    stub:
    "touch ${id}.vcfstats.xls"
}