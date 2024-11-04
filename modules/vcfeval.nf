process vcfeval {
    
    cpus params.CPU1
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(lib)
    tuple val(id), path(vcf) //demo.pf.bwa.gatk.vcf.gz

    output:
    tuple val(id), path("*.xls")
    tag "$id"
    publishDir "${params.outdir}/report/$id", saveAs: {"${id}.${lib}.evaluation.xls"}, mode: 'link'
    publishDir "${params.outdir}/report/$id", saveAs: {"${id}.evaluation.xls"}, mode: 'link', enabled: params.stLFR_only
    // publishDir (
    //   path: "${params.outdir}/report/", 
    //   saveAs: { fn ->
    //     if (fn.contains("lariat") && fn.contains("gatk")) { "${id}.52.lariat.gatk.vcf_eval.xls" }
    //     else if (fn.contains("lariat") && fn.contains("dv")) { "${id}.52.lariat.dv.vcf_eval.xls" }
    //     else if (fn.contains("bwa") && fn.contains("gatk")) { "${id}.52.bwa.gatk.vcf_eval.xls" }
    //     else { "${id}.52.bwa.dv.vcf_eval.xls" }
    // }
    // )
    
    script:
    def vcf = vcf.first()
    def prefix = "${vcf.getBaseName(2)}.${lib}"
    def benchmark = "${params.DB}/${params.ref}/${params.ref}.${params.std}.vcf.gz"
    def bed = "${params.DB}/${params.ref}/${params.ref}.${params.std}.bed"
    def sdf = "${params.DB}/${params.ref}/${params.ref}.SDF"
    """
    ${params.BIN}bcftools view -O z --type snps $vcf > ${id}.snp.vcf.gz
    ${params.BIN}tabix -p vcf -f ${id}.snp.vcf.gz

    ${params.BIN}bcftools view -O z --type indels $vcf > ${id}.indel.vcf.gz
    ${params.BIN}tabix -p vcf -f ${id}.indel.vcf.gz

    ${params.BIN}bcftools view -O z --type snps $benchmark > bench.snp.vcf.gz
    ${params.BIN}tabix -p vcf -f bench.snp.vcf.gz

    ${params.BIN}bcftools view -O z --type indels $benchmark > bench.indel.vcf.gz
    ${params.BIN}tabix -p vcf -f bench.indel.vcf.gz

    rm -rf ${id}_*_GIAB_v4

    ###vcf_eval
    ${params.BIN}rtg vcfeval \\
    -b bench.snp.vcf.gz \\
    -c ${id}.snp.vcf.gz \\
    -e $bed \\
    -t $sdf \\
    -T ${task.cpus} \\
    -o ${id}_snp_GIAB_v4

    ${params.BIN}rtg vcfeval \\
    -b bench.indel.vcf.gz \\
    -c ${id}.indel.vcf.gz \\
    -e $bed \\
    -t $sdf \\
    -T ${task.cpus} \\
    -o ${id}_indel_GIAB_v4

    head -1 ${id}_snp_GIAB_v4/summary.txt \\
      | awk 'BEGIN{OFS="\\t"}{print \$1,\$3,\$4,\$5,\$6,\$7,\$8}' \\
      > ${prefix}.evaluation.xls

    sed -i "s/None/SNP/g" ${id}_snp_GIAB_v4/summary.txt
    cat ${id}_snp_GIAB_v4/summary.txt \\
      | tail -1 \\
      | awk 'BEGIN{OFS="\\t"}{print \$1,\$3,\$4,\$5,\$6,\$7,\$8}' \\
      >> ${prefix}.evaluation.xls

    sed -i "s/None/Indel/g" ${id}_indel_GIAB_v4/summary.txt
    cat ${id}_indel_GIAB_v4/summary.txt \\
      | tail -1 \\
      | awk 'BEGIN{OFS="\\t"}{print \$1,\$3,\$4,\$5,\$6,\$7,\$8}' \\
      >> ${prefix}.evaluation.xls 
    rm bench.*.vcf.gz   
    """
    stub:
    "touch ${id}.${lib}.evaluation.xls"
}

process eachstat_vcf {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(vcf) //demo.pf.bwa.gatk.vcf.gz

    output:
    tuple val(id), path("${id}.varianttable.xls")

    publishDir "${params.outdir}/report/$id", mode: 'link'
    // publishDir (
    //   path: "${params.outdir}/report/", 
    //   saveAs: { fn ->
    //     if (fn.contains("lariat") && fn.contains("gatk")) { "${id}.52.lariat.gatk.vcf_eval.xls" }
    //     else if (fn.contains("lariat") && fn.contains("dv")) { "${id}.52.lariat.dv.vcf_eval.xls" }
    //     else if (fn.contains("bwa") && fn.contains("gatk")) { "${id}.52.bwa.gatk.vcf_eval.xls" }
    //     else { "${id}.52.bwa.dv.vcf_eval.xls" }
    // }
    // )
    
    script:
    def vcf = vcf.first()
    """
    ${params.BIN}bcftools view --types snps $vcf > snps.vcf
    ${params.BIN}bcftools view --types snps $vcf > indels.vcf
    ${params.BIN}python ${params.SCRIPT}/stat/eachstat_vcf.py $id > ${id}.varianttable.xls 
    rm *vcf
    """
}
process variant_fix {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(vcf), path(sv), path(cnv) //demo.pf.bwa.gatk.vcf.gz

    output:
    path "${id}.var.xls"

    tag "$id"
    publishDir "${params.outdir}/report/$id/", mode: 'link', saveAs: {"${id}.varianttable.xls"}
    // publishDir (
    //   path: "${params.outdir}/report/", 
    //   saveAs: { fn ->
    //     if (fn.contains("lariat") && fn.contains("gatk")) { "${id}.52.lariat.gatk.vcf_eval.xls" }
    //     else if (fn.contains("lariat") && fn.contains("dv")) { "${id}.52.lariat.dv.vcf_eval.xls" }
    //     else if (fn.contains("bwa") && fn.contains("gatk")) { "${id}.52.bwa.gatk.vcf_eval.xls" }
    //     else { "${id}.52.bwa.dv.vcf_eval.xls" }
    // }
    // )
    
    script:
    def vcf = vcf.first()
    """
    ${params.BIN}perl ${params.SCRIPT}/stat/variant_fix.pl $vcf $sv $cnv $id
    """
}