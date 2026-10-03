workflow WF_callvariants {
  take:
  ch_bam
  
  main:
  def ch_lariat = params.align_tool
  
  if (params.var_tool.contains("dv")) {
    if (params.use_megabolt) {dvMegabolt(ch_lariat, ch_bam).set {ch_mergevcf}}
    else {
      inferDvSex(ch_bam).set {ch_dv_sex}
      deepvariant(ch_lariat, ch_bam.join(ch_dv_sex)).set {ch_mergevcf}
    }
    
  } else if (params.var_tool.contains("gatk")) {
    if (params.use_megabolt) {
        if (params.run_vqsr) {
            vqsrMegabolt(ch_lariat, hcMegabolt(ch_lariat, ch_bam)).set {ch_mergevcf}
        } else {
            hcMegabolt(ch_lariat, ch_bam).set {ch_mergevcf}
        }
    } else {
        if (params.split_by_intervals) {
            gatk_interval().splitText().map { it.trim() }.collect().set {intervals}
            hcSplit(ch_lariat, ch_bam, intervals).set {ch_mergevcfSplit}
            gatherVcfsHc(ch_lariat, ch_mergevcfSplit.groupTuple()).set {ch_mergevcf}
            if (params.run_vqsr) {
                vqsrSnp(ch_lariat, ch_mergevcf).set {ch_vqsrsnp}
                vqsrIndel(ch_lariat, ch_mergevcf).set {ch_vqsrindel}
                gatherVcfsVqsr(ch_lariat, ch_vqsrsnp.join(ch_vqsrindel)).set{ch_mergevcf}
            }
            
        } else {
            hc(ch_lariat, ch_bam).set {ch_mergevcf}
            if (params.run_vqsr) {
                vqsrSnp(ch_lariat, hc(ch_lariat, ch_bam)).set {ch_vqsrsnp}
                vqsrIndel(ch_lariat, hc(ch_lariat, ch_bam)).set {ch_vqsrindel}
                gatherVcfsVqsr(ch_lariat, ch_vqsrsnp.join(ch_vqsrindel)).set{ch_mergevcf}
            } 
        }
    }
  }


  emit:
  ch_mergevcf
}
//megabolt
process hcMegabolt {
    label 'megabolt'
    cpus params.cpu3
    memory params.MEM2 + "g"
    
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

    ${params.MEGABOLT_RUNIT}  -l\$(basename \$(dirname \$PWD))_\$(basename \$PWD).${task.process}.${task.index} ${params.MEGABOLT} \\
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

    ${params.MEGABOLT_RUNIT}  -l\$(basename \$(dirname \$PWD))_\$(basename \$PWD).${task.process}.${task.index} ${params.MEGABOLT} \\
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

      gatk_tmp=/tmp/gatk_snp_${id}_\${BASHPID}
      trap 'rm -rf \${gatk_tmp}' EXIT
      mkdir -p \${gatk_tmp}
      ${params.BIN}gatk --java-options "-Xmx${task.memory.giga}g" \\
        VariantRecalibrator \\
        -V raw.snp.vcf.gz --tmp-dir \${gatk_tmp} \\
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

      rm -f snp.vcf.gz raw.snp.vcf.gz*
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

      gatk_tmp=/tmp/gatk_indel_${id}_\${BASHPID}
      trap 'rm -rf \${gatk_tmp}' EXIT
      mkdir -p \${gatk_tmp}
      ${params.BIN}gatk --java-options "-Xmx${task.memory.giga}g" \\
        VariantRecalibrator \\
        -V raw.indel.vcf.gz --tmp-dir \${gatk_tmp} \\
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

      rm -f indel.vcf.gz raw.indel.vcf.gz*
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
  cpus params.CPU0
  memory params.MEM0 + "g"

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
    
    input:
    val(aligner)
    tuple val(id), path(bam)

    output:
    tuple val(id), path("${id}.${aligner}.dv.vcf.gz*") //demo.lariat.dv.vcf.gz

    tag "$id"
    publishDir "${params.outdir}/$id/align/", mode: 'link'
 
    script:
    def bam = bam.first()
    def ref = params.ref.startsWith('/') ? params.ref : "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    """
    ${params.MEGABOLT_EXPORT}

    ${params.MEGABOLT_RUNIT}  -l\$(basename \$(dirname \$PWD))_\$(basename \$PWD).${task.process}.${task.index} ${params.MEGABOLT} \\
      --type haplotypecaller --haplotypecaller-input $bam --deepvariant 1 --fast-model 0 --ref $ref \\
      --outputdir .

    mv output/output.dv.vcf.gz ${id}.${aligner}.dv.vcf.gz
    mv output/output.dv.vcf.gz.tbi ${id}.${aligner}.dv.vcf.gz.tbi
    """
}
process inferDvSex {
    cpus 1
    memory params.MEM0 + "g"

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path("${id}.dv_sex.tsv")

    tag "$id"
    publishDir "${params.outdir}/$id/align/", mode: 'link'

    script:
    def input_bam = bam.find { it.toString().endsWith('.bam') } ?: bam.first()
    """
    case "${params.dv_sex_mode}" in
      auto)
        chr_y_length=`${params.BIN}samtools idxstats $input_bam | awk -v contig="${params.dv_sex_chrY_contig}" '\$1 == contig {print \$2; exit}'`
        if [ -z "\${chr_y_length}" ] || [ "\${chr_y_length}" -eq 0 ]; then
          echo "Cannot determine sex for $id: ${params.dv_sex_chrY_contig} is absent from $input_bam" >&2
          exit 2
        fi
        chr_y_bases=`${params.BIN}samtools depth -r "${params.dv_sex_chrY_contig}" $input_bam | awk '{sum += \$3} END {print sum + 0}'`
        chr_y_mean_depth=`awk -v bases="\${chr_y_bases}" -v contig_length="\${chr_y_length}" 'BEGIN {printf "%.6f", bases / contig_length}'`
        if awk -v depth="\${chr_y_mean_depth}" -v threshold="${params.dv_female_max_chrY_mean_depth}" 'BEGIN {exit !(depth < threshold)}'; then
          sex=female
        else
          sex=male
        fi
        ;;
      male|female)
        sex="${params.dv_sex_mode}"
        chr_y_mean_depth=NA
        ;;
      *)
        echo "Invalid dv_sex_mode: ${params.dv_sex_mode} (expected auto, male, or female)" >&2
        exit 2
        ;;
    esac
    printf 'sex\\tchrY_mean_depth\\tthreshold\\tmode\\n%s\\t%s\\t%s\\t%s\\n' "\${sex}" "\${chr_y_mean_depth}" "${params.dv_female_max_chrY_mean_depth}" "${params.dv_sex_mode}" > ${id}.dv_sex.tsv
    """
    stub:
    "printf 'sex\\tchrY_mean_depth\\tthreshold\\tmode\\nfemale\\t0.000000\\t${params.dv_female_max_chrY_mean_depth}\\tstub\\n' > ${id}.dv_sex.tsv"
}

process deepvariant {
    cpus params.cpu3
    memory params.deepvariant_memory
    maxForks 1
    container "${params.dv_container}"
    containerOptions "--env PATH=/opt/deepvariant/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/sbin:/bin:/usr/bin"

    input:
    val(aligner)
    tuple val(id), path(bam), path(sex)

    output:
    tuple val(id), path("${id}.*.dv.vcf.gz*") //demo.lariat.dv.vcf.gz

    tag "$id"
    publishDir "${params.outdir}/$id/align/", mode: 'link'
    beforeScript "export PATH=/opt/deepvariant/bin:\$PATH"
 
    script:
    def bam = bam.find { it.toString().endsWith('.bam') } ?: bam.first()
    def sex_file = sex
    def ref = params.ref.startsWith('/') ? params.ref : "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    def pangenome = params.dv_pangenome ?: "${params.DB}/hg38/panGenome/hprc-v1.1-mc-grch38.gbz"
    def ver = "dv"
    def outvcf = bam.toString().contains("pf") ? "${id}.pf.bwa.${ver}.vcf.gz" : "${id}.${aligner}.${ver}.vcf.gz"
    def gbz_shm = params.dv_gbz_shm_size_gb ? "--gbz_shared_memory_size_gb ${params.dv_gbz_shm_size_gb}" : ""
    def male_haploid_args = (params.dv_haploid_contigs && params.dv_haploid_contigs != "") ?
        "--haploid_contigs=${params.dv_haploid_contigs}" : ""
    def default_par_bed = pangenome.contains("/") ?
        "${pangenome.substring(0, pangenome.lastIndexOf('/'))}/GRCh38_PAR.bed" : ""
    def par_bed = (params.dv_par_regions_bed && params.dv_par_regions_bed != "") ?
        params.dv_par_regions_bed : default_par_bed
    def par_args = par_bed ? "--par_regions_bed=${par_bed}" : ""
    def base_make_args = params.dv_make_examples_extra_args ?: ""
    def male_make_args = [base_make_args, male_haploid_args, par_args].findAll { it }.join(",")
    def make_examples_flag = ""
    def postprocess_flag = (params.dv_postprocess_variants_extra_args && params.dv_postprocess_variants_extra_args != "") ?
        "--postprocess_variants_extra_args '${params.dv_postprocess_variants_extra_args}'" : ""
    def customized_model_flag = (params.dv_customized_model && params.dv_customized_model != "") ?
        "--customized_model=${params.dv_customized_model}" : ""
    def optional_args = [gbz_shm, customized_model_flag, make_examples_flag, postprocess_flag].findAll { it }.join(" \\\n      ")
    """
    dv_tmp=/tmp/dv_${id}_\${BASHPID}
    export TEST_TMPDIR=\${dv_tmp}/bazel
    export HOME=\${dv_tmp}/home
    trap 'rm -rf \${dv_tmp}' EXIT
    rm -rf \${dv_tmp}
    mkdir -p \${TEST_TMPDIR} \${HOME} \${dv_tmp}/intermediate

    sex=`awk 'NR == 2 {print \$1; exit}' $sex_file`
    case "\${sex}" in
      male)
        make_examples_args='${male_make_args}'
        ;;
      female)
        make_examples_args='${base_make_args}'
        ;;
      *)
        echo "Invalid sex value in $sex_file: \${sex}" >&2
        exit 2
        ;;
    esac
    make_examples_flag=()
    if [ -n "\${make_examples_args}" ]; then
      make_examples_flag=(--make_examples_extra_args "\${make_examples_args}")
    fi

    ${params.dv_binary_path} \\
      --model_type WGS \\
      --ref $ref \\
      --reads $bam \\
      --pangenome $pangenome \\
      --output_vcf $outvcf \\
      --num_shards ${task.cpus} \\
      --intermediate_results_dir \${dv_tmp}/intermediate \\
      $optional_args \\
      \${make_examples_flag[@]}
    """
    stub:
    def ver = "dv"
    def outvcf = bam.toString().contains("pf") ? "${id}.pf.bwa.${ver}.vcf.gz" : "${id}.${aligner}.${ver}.vcf.gz"
    "touch $outvcf"
}
