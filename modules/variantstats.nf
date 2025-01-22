process variant_stats {
    
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("*.xls"), emit: varianttable

    tag "$id"
    // publishDir (
    //   path: "${params.outdir}/report/", 
    //   saveAs: { fn ->
    //     if (fn.contains("lariat") && fn.contains("gatk")) { "${id}.51.lariat.gatk.variant_stats.xls" }
    //     else if (fn.contains("lariat") && fn.contains("dv")) { "${id}.51.lariat.dv.variant_stats.xls" }
    //     else if (fn.contains("bwa") && fn.contains("gatk")) { "${id}.51.bwa.gatk.variant_stats.xls" }
    //     else { "${id}.51.bwa.dv.variant_stats.xls" }
    // }
    // )

    script:
    def vcffile = vcf.first()
    def prefix = "${vcffile.getBaseName(2)}"
    """
    # <vcf> <out> <bcftools> <sample>
    perl ${params.SCRIPT}/stat/eachstat_vcf.pl $vcffile ${prefix}.variant_stats.xls ${params.BIN}bcftools $prefix
    """
}