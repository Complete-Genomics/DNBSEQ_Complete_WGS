
//megabolt
process hcMegabolt {
    label 'megabolt'
    cpus params.cpu3
    memory params.MEM2 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.boltq} ${params.extraCluOpt}"
    
    input:
    val(aligner)
    tuple val(id), path(bam) //demo.pf.bwa.merge.bam

    output:
    tuple val(id), path("${id}.${aligner}.megaboltHc?.vcf.gz*")

    tag "$id"
    // publishDir "${params.outdir}/$id/align/"
 
    script:
    def bam = bam.first()
    def gatk = params.gatk_version == "v4" ? "--hc4 1 --gatk4 1" : ""
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    def outprefix = params.gatk_version == "v4" ? "${id}.${aligner}.megaboltHc4" : "${id}.${aligner}.megaboltHc3"
    """
    ${params.MEGABOLT_EXPORT}
    
    hapmap=`ls ${params.DB}/${params.ref}/gatk/*hapmap*.vcf.gz`
    dbsnp=`ls ${params.DB}/${params.ref}/gatk/*dbsnp*.vcf.gz`
    mills=`ls ${params.DB}/${params.ref}/gatk/Mills*.vcf.gz`
    kgsnp=`ls ${params.DB}/${params.ref}/gatk/1000G*snps*.vcf.gz`

    ${params.MEGABOLT_RUNIT} -l${task.process}.${task.index} ${params.MEGABOLT}                           \\
      --type haplotypecaller $gatk \\
      --haplotypecaller-input $bam \\
      --stand-call-conf 10 --ref $ref \\
      --vcf \$dbsnp \\
      --knownSites \$dbsnp \\
      --knownSites \$mills \\
      --knownSites \$kgsnp \\
      --outputdir . 

    mv output/output.hc*.vcf.gz ${outprefix}.vcf.gz
    mv output/output.hc*.vcf.gz.tbi ${outprefix}.vcf.gz.tbi
    """
    stub:
    "touch ${id}.${aligner}.megaboltHc?.vcf.gz"
}
process vqsrMegabolt {
    label 'megabolt'
    cpus params.cpu3
    memory params.MEM2 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.boltq} ${params.extraCluOpt}"
    
    input:
    val(aligner)
    tuple val(id), path(vcf) //demo.pf.megaboltHc.vcf.gz

    output:
    tuple val(id), path("${id}.${aligner}.gatk?.vcf.gz*")

    tag "$id"
    publishDir "${params.outdir}/$id/align/", pattern: "*vcf.gz", saveAs: {"${id}.${aligner}.megaboltVqsr*.vcf.gz"}, mode: 'link'
    publishDir "${params.outdir}/$id/align/", pattern: "*tbi", saveAs: {"${id}.${aligner}.megaboltVqsr*.vcf.gz.tbi"}, mode: 'link'
 
    script:
    def vcf = vcf.first()
    def gatk = params.gatk_version == "v4" ? "--hc4 1 --gatk4 1" : ""
    def outprefix = params.gatk_version == "v4" ? "${id}.${aligner}.gatk4" : "${id}.${aligner}.gatk3"
    """
    ${params.MEGABOLT_EXPORT}
    
    hapmap=`ls ${params.DB}/${params.ref}/gatk/*hapmap*.vcf.gz`
    dbsnp=`ls ${params.DB}/${params.ref}/gatk/*dbsnp*.vcf.gz`
    mills=`ls ${params.DB}/${params.ref}/gatk/Mills*.vcf.gz`
    kgsnp=`ls ${params.DB}/${params.ref}/gatk/1000G*snps*.vcf.gz`
    omni=`ls ${params.DB}/${params.ref}/gatk/*omni*.vcf.gz`

    ${params.MEGABOLT_RUNIT} -l${task.process}.${task.index} ${params.MEGABOLT}                           \\
      --type vqsr --vqsr-input $vcf $gatk \\
      --resource-hapmap \$hapmap \\
      --resource-omni \$omni \\
      --resource-1000G \$kgsnp \\
      --resource-dbsnp \$dbsnp \\
      --resource-mills \$mills \\
      --outputdir . 

    mv output/output.vqsr.vcf.gz ${outprefix}.vcf.gz
    mv output/output.vqsr.vcf.gz.tbi ${outprefix}.vcf.gz.tbi
    """
    stub:
    "touch ${id}.${aligner}.gatk?.vcf.gz"
}
//no megabolt
process hc {
    cpus params.cpu3
    memory params.MEM2 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(aligner)
    tuple val(id), path(bam) //demo.pf.bwa.merge.bam

    output:
    tuple val(id), path("${id}.${aligner}.*.vcf.gz*")

    tag "$id"
    publishDir "${params.outdir}/$id/align/"
 
    script:
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    bam = bam.first()
    cmd = """
    hapmap=`ls ${params.DB}/${params.ref}/gatk/*hapmap*.vcf.gz`
    dbsnp=`ls ${params.DB}/${params.ref}/gatk/*dbsnp*.vcf.gz`
    mills=`ls ${params.DB}/${params.ref}/gatk/Mills*.vcf.gz`
    kgsnp=`ls ${params.DB}/${params.ref}/gatk/1000G*snps*.vcf.gz`
    kgindel=`ls ${params.DB}/${params.ref}/gatk/1000G*indels*.vcf.gz`
    omni=`ls ${params.DB}/${params.ref}/gatk/*omni*.vcf.gz`
    """
    if (params.gatk_version == "v4") {
      cmd += """
      ${params.BIN}gatk --java-options "-Xmx${task.memory.giga}g" \\
        HaplotypeCaller \\
        -R $ref \\
        -I $bam \\
        --dbsnp \$dbsnp \\
        -O ${id}.${aligner}.hc4.vcf.gz
      ${params.BIN}tabix ${id}.${aligner}.hc4.vcf.gz
      """
    } else if (params.gatk_version == "v3") {
      cmd += """
      ${params.BIN}gatk3 -Xmx${task.memory.giga}g -Djava.io.tmpdir=tmpdir \\
        -T HaplotypeCaller \\
        -nct ${task.cpus}  \\
        -R $ref \\
        -I $bam \\
        --dbsnp \$dbsnp \\
        -o ${id}.${aligner}.hc3.vcf.gz
      ${params.BIN}tabix ${id}.${aligner}.hc3.vcf.gz
      """
    }
    return cmd
    stub:
    "touch ${id}.${aligner}.*.vcf.gz"
}
process vqsrSnp {
    cpus params.cpu3
    memory params.MEM2 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(aligner)
    tuple val(id), path(vcf) //id.aligner.hc*.vcf.gz

    output:
    tuple val(id), path("${id}.${aligner}.vqsr?.snp.vcf.gz*")

    tag "$id"
    publishDir "${params.outdir}/$id/align/"
 
    script:
    vcf = vcf.first()
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    cmd = """
    hapmap=`ls ${params.DB}/${params.ref}/gatk/*hapmap*.vcf.gz`
    dbsnp=`ls ${params.DB}/${params.ref}/gatk/*dbsnp*.vcf.gz`
    mills=`ls ${params.DB}/${params.ref}/gatk/Mills*.vcf.gz`
    kgsnp=`ls ${params.DB}/${params.ref}/gatk/1000G*snps*.vcf.gz`
    kgindel=`ls ${params.DB}/${params.ref}/gatk/1000G*indels*.vcf.gz`
    omni=`ls ${params.DB}/${params.ref}/gatk/*omni*.vcf.gz`
    """
    if (params.gatk_version == "v4") {
      cmd += """
      ${params.BIN}gatk --java-options "-Xmx${task.memory.giga}g" \\
        SelectVariants -R $ref -V $vcf -select-type SNP --exclude-non-variants -O raw.snp.vcf.gz

      mkdir tmpdir
      ${params.BIN}gatk --java-options "-Xmx${task.memory.giga}g" \\
        VariantRecalibrator \\
        -V raw.snp.vcf.gz --tmp-dir tmpdir \\
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 \$hapmap \\
        --resource:omni,known=false,training=true,truth=true,prior=12.0 \$omni \\
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 \$kgsnp \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \$dbsnp \\
        -an DP -an QD -an FS -an SOR -an ReadPosRankSum \\
        -mode SNP \\
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\
        -O recalibrate.recal \\
        --tranches-file recalibrate.tranches

      ${params.BIN}gatk --java-options "-Xmx${task.memory.giga}g" \\
        ApplyVQSR \\
        -V raw.snp.vcf.gz \\
        --recal-file recalibrate.recal \\
        -mode SNP \\
        -O snp.vcf.gz

      ${params.BIN}gatk --java-options "-Xmx${task.memory.giga}g" \\
        SelectVariants -R $ref -V snp.vcf.gz --exclude-filtered -O ${id}.${aligner}.vqsr4.snp.vcf.gz

      rm -rf tmpdir snp.vcf.gz raw.snp.vcf.gz*
      """
    } else if (params.gatk_version == "v3") {
      cmd += """
      ${params.BIN}gatk3 -Xmx${task.memory.giga}g -Djava.io.tmpdir=tmpdir \\
        -T SelectVariants \\
        -R $ref -V $vcf -selectType SNP --excludeNonVariants -o raw.snp.vcf.gz

      ${params.BIN}gatk3 -Xmx${task.memory.giga}g -Djava.io.tmpdir=tmpdir \\
        -T  VariantRecalibrator \\
        -R $ref \\
        -nt ${task.cpus} \\
        -input raw.snp.vcf.gz \\
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 \$hapmap \\
        -resource:omni,known=false,training=true,truth=true,prior=12.0 \$omni \\
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 \$kgsnp \\
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \$dbsnp \\
        -an DP -an QD -an FS -an SOR -an ReadPosRankSum \\
        -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\
        -recalFile recalibrate.recal -tranchesFile recalibrate.tranches -rscriptFile recalibrate_plots.R

      ${params.BIN}gatk3 -Xmx${task.memory.giga}g -Djava.io.tmpdir=tmpdir \\
        -T ApplyRecalibration \\
        -R $ref \\
        -input raw.snp.vcf.gz \\
        -mode SNP \\
        --ts_filter_level 99.9 -recalFile recalibrate.recal -tranchesFile recalibrate.tranches \\
        -o snp.vcf.gz   

      ${params.BIN}gatk3 -Xmx${task.memory.giga}g -Djava.io.tmpdir=tmpdir \\
        -T SelectVariants -R $ref -V snp.vcf.gz --excludeFiltered -o ${id}.${aligner}.vqsr3.snp.vcf.gz \\

      rm raw.snp.vcf.gz* snp.vcf.gz
      """
    }
    return cmd
    stub:
    "touch ${id}.${aligner}.vqsr?.snp.vcf.gz"
}
process vqsrIndel {
    cpus params.cpu3
    memory params.MEM2 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(aligner)
    tuple val(id), path(vcf) //demo.pf.bwa.merge.bam

    output:
    tuple val(id), path("${id}.${aligner}.vqsr?.indel.vcf.gz*")

    tag "$id"
    publishDir "${params.outdir}/$id/align/"
 
    script:
    vcf = vcf.first()
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    cmd = """
    hapmap=`ls ${params.DB}/${params.ref}/gatk/*hapmap*.vcf.gz`
    dbsnp=`ls ${params.DB}/${params.ref}/gatk/*dbsnp*.vcf.gz`
    mills=`ls ${params.DB}/${params.ref}/gatk/Mills*.vcf.gz`
    kgsnp=`ls ${params.DB}/${params.ref}/gatk/1000G*snps*.vcf.gz`
    kgindel=`ls ${params.DB}/${params.ref}/gatk/1000G*indels*.vcf.gz`
    omni=`ls ${params.DB}/${params.ref}/gatk/*omni*.vcf.gz`
    """
    if (params.gatk_version == "v4") {
      cmd += """
      ${params.BIN}gatk --java-options "-Xmx${task.memory.giga}g" \\
        SelectVariants -R $ref -V $vcf -select-type INDEL --exclude-non-variants -O raw.indel.vcf.gz

      mkdir tmpdir
      ${params.BIN}gatk --java-options "-Xmx${task.memory.giga}g" \\
        VariantRecalibrator \\
        -V raw.indel.vcf.gz --tmp-dir tmpdir \\
        -resource:mills,known=true,training=true,truth=true,prior=12.0 \$mills \\
        -an DP -an QD -an FS -an SOR -an ReadPosRankSum \\
        -mode INDEL \\
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\
        -O recalibrate.recal \\
        --tranches-file recalibrate.tranches

      ${params.BIN}gatk --java-options "-Xmx${task.memory.giga}g" \\
        ApplyVQSR \\
        -V raw.indel.vcf.gz \\
        --recal-file recalibrate.recal \\
        -mode INDEL \\
        -O indel.vcf.gz

      ${params.BIN}gatk --java-options "-Xmx${task.memory.giga}g" \\
        SelectVariants -R $ref -V indel.vcf.gz --exclude-filtered -O ${id}.${aligner}.vqsr4.indel.vcf.gz

      rm -rf tmpdir indel.vcf.gz raw.indel.vcf.gz* 
      """
    } else if (params.gatk_version == "v3") {
      cmd += """
      ${params.BIN}gatk3 -Xmx${task.memory.giga}g -Djava.io.tmpdir=tmpdir \\
        -T SelectVariants \\
        -R $ref -V $vcf -selectType INDEL --excludeNonVariants -o raw.indel.vcf.gz

      ${params.BIN}gatk3 -Xmx${task.memory.giga}g -Djava.io.tmpdir=tmpdir \\
        -T  VariantRecalibrator \\
        -R $ref \\
        -nt ${task.cpus} \\
        -input raw.indel.vcf.gz \\
        -resource:mills,known=true,training=true,truth=true,prior=12.0 \$mills \\
        -an DP -an QD -an FS -an SOR -an ReadPosRankSum \\
        -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 \\
        -recalFile recalibrate.recal -tranchesFile recalibrate.tranches -rscriptFile recalibrate_plots.R

      ${params.BIN}gatk3 -Xmx${task.memory.giga}g -Djava.io.tmpdir=tmpdir \\
        -T ApplyRecalibration \\
        -R $ref \\
        -input raw.indel.vcf.gz \\
        -mode INDEL \\
        --ts_filter_level 99.9 -recalFile recalibrate.recal -tranchesFile recalibrate.tranches \\
        -o indel.vcf.gz  
      
      ${params.BIN}gatk3 -Xmx${task.memory.giga}g -Djava.io.tmpdir=tmpdir \\
        -T SelectVariants -R $ref -V indel.vcf.gz --excludeFiltered -o ${id}.${aligner}.vqsr3.indel.vcf.gz \\

      rm raw.indel.vcf.gz* indel.vcf.gz*
      """
    }
    return cmd
    stub:
    "touch ${id}.${aligner}.vqsr?.indel.vcf.gz"
}
//run hc split
process gatk_interval {
  executor = 'local'
	container false

  cpus params.CPU0
  memory params.MEM0 + "g"
  clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

  output:
  file("txt")

  script:
  def fai = "${params.DB}/${params.ref}/reference/${params.ref}.fa.fai"
  """
  #!/usr/bin/env python

  import os
  chrs = ["chr" + str(i) for i in range(1,23)]
  chrs.append("chrX")

  f = open("$fai")
  g2 = open("chrOthers.bed", 'w')
  for line in f:
    chr, l, _, _, _ = line.rstrip().split()
    if chr in chrs:
      g = open(chr + ".bed", 'w')
      g.write(chr + "\\t0\\t" + l + "\\n")
      g.close()
    else:
      g2.write(chr + "\\t0\\t" + l + "\\n")
  f.close()
  g2.close()

  g3 = open("txt", 'w')
  for bed in os.listdir("."):
    if not bed.endswith("bed"):
      continue
    g3.write(os.path.abspath(bed) + "\\n")
  g3.close()

  #head -22 $fai | awk '{print \$1}' > txt
  #tail -n +23 $fai |awk '{printf "%s -L ", \$1}' |sed 's/-L \$//' >> txt
  """
  stub:
  "touch txt"
}
process hcSplit {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(aligner)
    tuple val(id), path(bam) //demo.pf.bwa.merge.bam
    each bed

    output:
    tuple val(id), path("${id}.${aligner}.*.vcf.gz")

    tag "$id, ${file(bed).getBaseName()}"
 
    script:
    bam = bam.first()
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    def chr1 = file(bed).getBaseName()
    cmd = """
    hapmap=`ls ${params.DB}/${params.ref}/gatk/*hapmap*.vcf.gz`
    dbsnp=`ls ${params.DB}/${params.ref}/gatk/*dbsnp*.vcf.gz`
    mills=`ls ${params.DB}/${params.ref}/gatk/Mills*.vcf.gz`
    kgsnp=`ls ${params.DB}/${params.ref}/gatk/1000G*snps*.vcf.gz`
    kgindel=`ls ${params.DB}/${params.ref}/gatk/1000G*indels*.vcf.gz`
    omni=`ls ${params.DB}/${params.ref}/gatk/*omni*.vcf.gz`
    """
    if (params.gatk_version == "v4") {
      cmd += """
      ${params.BIN}gatk --java-options "-Xmx${task.memory.giga}g" \\
        HaplotypeCaller \\
        -R $ref \\
        -I $bam \\
        --dbsnp \$dbsnp \\
        -L $bed \\
        -O ${id}.${aligner}.${chr1}.hc4.vcf.gz
      """
    } else if (params.gatk_version == "v3") {
      cmd += """
      ${params.BIN}gatk3 -Xmx${task.memory.giga}g -Djava.io.tmpdir=tmpdir \\
        -T HaplotypeCaller \\
        -nct ${task.cpus}  \\
        -R $ref \\
        -I $bam \\
        --dbsnp \$dbsnp \\
        -L $bed \\
        -o ${id}.${aligner}.${chr1}.hc3.vcf.gz
      """
    }
    return cmd
    stub:
    "touch ${id}.${aligner}.${file(bed).getBaseName()}.vcf.gz"
}
process gatherVcfsHc {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(aligner)
    tuple val(id), path(vcfs) //${id}.${aligner}.${chr1}.hc*.vcf.gz

    output:
    tuple val(id), path("${id}.${aligner}.*vcf.gz*")

    tag "$id"
    publishDir "${params.outdir}/$id/align/"
 
    script:
    vcf = vcfs.first()
    suffix = vcf.toString().tokenize('.')[-3] + ".vcf.gz"
    // println "$suffix"
    """
    vcfs=""
    for i in {1..22} X Others;do
      vcfs="\$vcfs -I ${id}.${aligner}.chr\$i.*.vcf.gz"
    done

    ${params.BIN}gatk \\
      GatherVcfs \\
      \$vcfs \\
      -O ${id}.${aligner}.${suffix}

    ${params.BIN}tabix ${id}.${aligner}.${suffix}
    """
    stub:
    "touch ${id}.${aligner}.${suffix}"
}
process gatherVcfsVqsr {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(aligner)
    tuple val(id), path(snpvcf), path(indelvcf) //${id}.${aligner}.vqsr3.indel.vcf.gz

    output:
    tuple val(id), path("${id}.${aligner}.*vcf.gz")

    tag "$id"
    publishDir "${params.outdir}/$id/align/"
 
    script:
    snpvcf = snpvcf.first()
    indelvcf = indelvcf.first()
    prefix = snpvcf.getBaseName(3)
    """
    ${params.BIN}bcftools concat -a $snpvcf $indelvcf -O b -o ${prefix}.vcf.gz
    """
    stub:
    "touch ${prefix}.vcf.gz"
}
process dvMegabolt {
  label 'megabolt'
    cpus params.cpu3
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.boltq} ${params.extraCluOpt}"
    
    input:
    val(aligner)
    tuple val(id), path(bam) //demo.stlfr.lariat.merge.bam

    output:
    tuple val(id), path("${id}.${aligner}.dv.vcf.gz*") //demo.lariat.dv.vcf.gz

    tag "$id"
    publishDir "${params.outdir}/$id/align/", mode: 'link'
 
    script:
    def bam = bam.first()
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    """
    dbsnp=`ls ${params.DB}/${params.ref}/gatk/*dbsnp*.vcf.gz`
    kgsnp=`ls ${params.DB}/${params.ref}/gatk/1000G*snps*.vcf.gz`
    mills=`ls ${params.DB}/${params.ref}/gatk/Mills*indel*.vcf.gz`
    
    ${params.MEGABOLT_EXPORT}

    ${params.MEGABOLT_RUNIT} -l${task.process}.${task.index} ${params.MEGABOLT}                           \\
      --type haplotypecaller --haplotypecaller-input $bam --deepvariant 1 --fast-model 0 --ref $ref \\
      --vcf \$dbsnp \\
      --knownSites \$kgsnp --knownSites \$mills \\
      --outputdir .

    mv output/output.dv.vcf.gz ${id}.${aligner}.dv.vcf.gz
    mv output/output.dv.vcf.gz.tbi ${id}.${aligner}.dv.vcf.gz.tbi
    """
    stub:
    "touch ${id}.${aligner}.dv.vcf.gz"
}
process deepvariantv16 {
    cpus params.cpu3
    memory params.MEM3 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(aligner)
    tuple val(id), path(bam) //demo.stlfr.lariat.merge.bam

    output:
    tuple val(id), path("${id}.*.dv.vcf.gz*") //demo.lariat.dv.vcf.gz

    tag "$id"
    publishDir "${params.outdir}/$id/align/", mode: 'link'
 
    script:
    def bam = bam.first()
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    def machine="${params.dv_machine}" // g400
    def model="${params.DB}/DV_model/dv1.6-mgi-${machine}.ckpt"
    def outvcf = bam.toString().contains("pf") ? "${id}.pf.bwa.dv.vcf.gz" : "${id}.${aligner}.dv.vcf.gz"
    """
    /opt/deepvariant/bin/run_deepvariant --model_type=WGS \\
      --ref=$ref \\
      --reads=$bam \\
      --output_vcf=$outvcf \\
      --logging_dir=log \\
      --customized_model=${model} \\
      --runtime_report=True \\
      --intermediate_results_dir="intermediate_results_dir" \\
      --num_shards=${task.cpus}\\
      --make_examples_extra_args="vsc_min_count_indels=1" 

    """
    stub:
    "touch ${id}.${aligner}.dv.vcf.gz"
}
