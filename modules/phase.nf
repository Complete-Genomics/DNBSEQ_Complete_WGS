workflow WF_phase {
    take:
    ch_input

    main:
    ch_input.map { id, bam, vcf -> [id, bam]}.set {ch_bam}
    ch_input.map { id, bam, vcf -> [id, vcf]}.set {ch_mergevcf}

    if (params.chr == 'all') {
        chrs = (1..22).collect { "chr$it" } + ["chrX", "chrY"]
    } else {
        chrs = [params.chr]
    }

    splitBam4phasing(ch_bam, chrs).set {ch_eachbam}
    splitvcf(ch_mergevcf, chrs).eachvcf.set {ch_eachvcf}

    vcfs = splitvcf.out.vcf.groupTuple()
    ch_eachbam.combine(ch_eachvcf, by: [0,1]).set {ch_eachchr}

    pvcfs = phase(ch_eachchr).phasedvcf.groupTuple()  
    lfs = phase.out.lf.groupTuple()  
    hbs = phase.out.hapblock.groupTuple()  
    stats = phase.out.stat.groupTuple()  

    if (params.ref == 'hg38' || params.ref.contains('GRCh38')) {
        phaseCat(vcfs.join(pvcfs).join(lfs).join(hbs).join(stats)).report.set {ch_phasereport}//report
        phaseCat.out.phasedvcf.set {ch_phasedvcf}
        phaseCat.out.hb.set {ch_hb}
    } else {
        getchrs().set { txt }
        chrs = txt.splitText().map { it.trim() }.collect()

        phaseCatRef(ch_lariat, ch_dv, txt, vcfs.join(pvcfs).join(lfs).join(hbs)).report.set {ch_phasereport}
        phaseCatRef.out.phasedvcf.set {ch_phasedvcf}
        phaseCatRef.out.hb.set {ch_hb}
    }

    emit:
    vcf = ch_phasedvcf
    hb = ch_hb
    report = ch_phasereport
}
process splitBam4phasing {
  
  cpus params.CPU0
  memory params.MEM0 + "g"
  clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
  

  input:
  tuple val(id), path(bam)
  each chr

  output:
  tuple val(id), val(chr), path("${id}.${params.align_tool}.${chr}.bam*"), emit: eachbam

  tag "$id"
  
  // publishDir "${params.outdir}/$id/align/alignsplit"

  script:
  def aligner = params.align_tool
  def bam = bam.first()
  """
  # ${params.BIN}samtools view -h -F 0x400 ${bam} ${chr} \\
  #  | awk -v OFS='\\t' '{if(\$1~/#/){split(\$1,a,"#"); if(a[2]!~/0_0_0/){sub(/BX:Z:.*/, "BX:Z:"a[2]); print \$0} }else{print}}' - \\
  #  | ${params.BIN}samtools view -bhS - > ${id}.${aligner}.${chr}.bam

  # ${params.BIN}samtools view -h -F 0x400 $bam $chr \\
  # | perl -ne '\$p1=index(\$_,"#");\$p2=index(\$_,"\\t",\$p1);if (\$p1>0 && \$p2>0) {\$p3=rindex(\$_,"\\tBX:Z:");if (\$p3>0){substr(\$_,\$p3)="\\tBX:Z:".substr(\$_,\$p1+1,\$p2-\$p1-1)."\\n"}}print' | ${params.BIN}samtools view -bhS - > ${id}.${aligner}.${chr}.bam
  
  
  ${params.BIN}samtools view -h -F 0x400 $bam $chr \\
  | perl -ne '\$p1=index(\$_,"#");\$p2=index(\$_,"\\t",\$p1);if (\$p1>=0 && \$p2>\$p1) {\$bx_tag="BX:Z:".substr(\$_,\$p1+1,\$p2-\$p1-1); \$p3=rindex(\$_,"\tBX:Z:"); if (\$p3>=0){substr(\$_,\$p3)="\\t\$bx_tag\\n"} else {chomp; \$_.="\\t\$bx_tag\\n"}} print' | ${params.BIN}samtools view -bhS - > ${id}.${aligner}.${chr}.bam

  ${params.BIN}samtools index ${id}.${aligner}.${chr}.bam 
  """
  stub:
  "touch ${id}.${aligner}.${chr}.bam "
}
process splitvcf {
  
  cpus params.CPU0
  memory params.MEM0 + "g"
  clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    
  input:
  tuple val(id), path(vcf)
  each chr

  output:
  tuple val(id), val(chr), path("*.vcf.gz"), emit: eachvcf
  tuple val(id), path("*.vcf.gz"), emit: vcf

  tag "$id"
  // publishDir "${params.outdir}/$id/align/alignsplit"

  script:
  def vcf = vcf.first()
  """
  zcat $vcf | awk -F'\\t' -v chr="$chr" '(\$0 ~ /^#/) || (\$1 == chr && \$10 ~ /^(1\\/1|0\\/1|1\\/2):/)' | ${params.BIN}bgzip > ${id}.${params.align_tool}.${params.var_tool}.${chr}.vcf.gz
  """
}

process getchrs {
    cpus params.CPU0
    memory params.MEM1 + "g"
    
    output:
    file("txt")

    script:
    def fai = "${params.ref}.fai"
    """
    awk '\$2 > ${params.cut_len} {print \$1}' $fai > txt 
    """
    stub:
    "touch txt"

}

process phase {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    

    input:
    tuple val(id), val(chr), path(bam), path(vcf) //demo.stlfr.gatk.chr1.vcf.gz

    output:
    tuple val(id), path("*.VCF.gz"), emit: phasedvcf
    tuple val(id), path("*.lf"), emit: lf
    tuple val(id), path("*.hapblock"), emit: hapblock
    tuple val(id), val(chr), path(bam), path("*.hapblock"), emit: svpre
    tuple val(id), path("*.hapcut_stat.txt"), emit: stat

    tag "$id, $chr"
    publishDir "${params.outdir}/$id/phase/phasesplit", mode: 'link'

    script:
    def bam = bam.first()
    def prefix = "${id}.${params.align_tool}.${params.var_tool}.${chr}"

    cmd = """
    #export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:${params.DB}/htslib
    ${params.BIN}gzip -dc $vcf > tmp.vcf
    ${params.BIN_HAPCUT2}extractHAIRS --indels 1 --10X 1 --bam $bam --VCF tmp.vcf --out ${prefix}.unlinked_frag # --triallelic 1

    python3 ${params.SCRIPT_HapCUT2}LinkFragments.py \\
        --bam $bam --VCF tmp.vcf --fragments ${prefix}.unlinked_frag --out ${prefix}.lf -d 100000

    ( ${params.BIN_HAPCUT2}HAPCUT2 --nf 1 --fragments ${prefix}.lf --vcf tmp.vcf --output ${prefix}.hapblock ) || \
    ( ${params.BIN_HAPCUT2}HAPCUT2 --nf 1 --fragments ${prefix}.lf --vcf tmp.vcf --output ${prefix}.hapblock --skip_prune 1 )

    ${params.BIN}bgzip ${prefix}.hapblock.phased.VCF
    rm tmp.vcf
    echo $chr > ${prefix}.hapcut_stat.txt
    """
    if (!params.ref.startsWith('/')) {
        def pv = "${params.DB}/${params.ref}/phasedvcf/${params.ref}.${params.std}.${chr}.vcf.gz"

        cmd += """
        if [ -s ${prefix}.hapblock ]; then
            python3 ${params.SCRIPT}/calculate_haplotype_statistics_CWX.py \\
            -h1 ${prefix}.hapblock -v1 ${prefix}.hapblock.phased.VCF.gz -v2 $pv --indels >> ${prefix}.hapcut_stat.txt
        fi
        """
    } else {
        cmd += """
        if [ -s ${prefix}.hapblock ]; then
            python3 ${params.SCRIPT}/calc_n50.py \\
            -h1 ${prefix}.hapblock -v1 ${prefix}.hapblock.phased.VCF.gz --indels >> ${prefix}.hapcut_stat.txt
        fi
        """
    }
        
    return cmd
}

process ideogram {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    when: params.ref == 'hg38' || params.ref.contains('GRCh38')

    input:
    tuple val(id), path(hapblock)

    output:
    tuple val(id), path("*png"), emit: png
    path("*svg")

    tag "$id"
    publishDir "${params.outdir}/report/$id/", mode: 'link'

    script:
    def fai = "${params.DB}/${params.ref}/reference/${params.ref}.fa.fai"
    """
    ${params.BIN}python ${params.SCRIPT}/band.py $hapblock
    /usr/local/miniconda3/envs/rideogram/bin/Rscript ${params.SCRIPT}/ideogram.R ${params.SCRIPT}/hg38.karyotype
    convert -crop 1250x1250+900+200 chromosome.png tmp
    mv tmp chromosome.png
    """
    stub:
    "touch chromosome.png chromosome.svg"
}
process cumuplot {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    
    when: params.ref == 'hg38' || params.ref.contains('GRCh38')

    input:
    tuple val(id), path(hapblock)

    output:
    tuple val(id), path("cumulative_coverage_plot.png")

    // tag "$id, $aligner, $varcaller"
    tag "$id"
    publishDir "${params.outdir}/report/$id/", mode: 'link'

    script:
    """
    /usr/local/miniconda3/envs/six/bin/python ${params.SCRIPT}/cumuplot.py $hapblock
    """
}

process hapKaryotype {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    

    input:
    // val(aligner)
    // val(varcaller)
    tuple val(id), path(hapblock)

    output:
    tuple val(id), path("*pdf"), emit: pdf
    path("*txt")

    // tag "$id, $aligner, $varcaller"
    tag "$id"
    publishDir "${params.outdir}/report/$id/", mode: 'link'

    script:
    def fai = "${params.DB}/${params.ref}/reference/${params.ref}.fa.fai"
    """
    # output: karyotype.id.genome.txt, karyotype.id.band.txt
    ${params.BIN}python ${params.SCRIPT}/stat/haplotype.karyotype.py $id $fai $hapblock

    ${params.BIN}Rscript ${params.SCRIPT}/stat/haplotype.karyotype.R karyotype.${id}.genome.txt karyotype.${id}.band.txt ${id}.haplotype.pdf .
    """
    stub:
    "touch ${id}.haplotype.pdf karyotype.${id}.band.txt"
}
process hapKaryotype_bak {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    
    // cache false
    input:
    // val(aligner)
    // val(varcaller)
    tuple val(id), path(hapblock)

    output:
    tuple val(id), path("*pdf")
    path("*txt")

    // tag "$id, $aligner, $varcaller"
    tag "$id"
    publishDir "${params.outdir}/report/$id/", pattern: "*pdf", saveAs: {"${id}.haplotype.perl.pdf"}, mode: 'link'
    publishDir "${params.outdir}/report/$id/", pattern: "karyotype.${id}.genome.txt", saveAs: {"karyotype.${id}.genome.perl.txt"}, mode: 'link'
    publishDir "${params.outdir}/report/$id/", pattern: "karyotype.${id}.band.txt", saveAs: {"karyotype.${id}.band.perl.txt"}, mode: 'link'

    script:
    def fa = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    """
    # output: karyotype.id.genome.txt, karyotype.id.band.txt
    ${params.BIN}perl ${params.SCRIPT}/stat/haplotype.karyotype.pl $id $hapblock . $fa .
    """
}
process circos {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    

    input:
    // val(aligner)
    // val(varcaller)
    tuple val(id), path(vcf), path(sv), path(cnv), path(bamfinal)

    output:
    path("${id}.circos.svg")
    path("${id}.circos.png")
    tuple val(id), path("*png"), emit: flg

    // cache false
    tag "$id"
    publishDir "${params.outdir}/report/$id/", mode: 'link'

    script:
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
	def vcf = vcf.first()
    // def prefix = "${id}.${aligner}.${varcaller}"
    //my ($name, $vcf, $fsv, $fcnv, $fbamfinal, $ref, $circos, $bcftools) = @ARGV;

    """    
    ${params.BIN}perl ${params.SCRIPT}/stat/circos.pl $id $vcf $sv $cnv $bamfinal $ref /usr/local/app/circos /usr/bin/bcftools
    ${params.BIN}circos -conf circos.${id}.conf -outputfile ${id}.circos.png -outputdir .

    """
}
process phaseCat {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    
    input:
    tuple val(id), path(vcfs), path(pvs), path(lfs), path(hbs), path(stats)

    output:
    tuple val(id), path("*hapblock"), emit: hb
    tuple val(id), path("*.hapcut_stat.txt"), emit: hapcutstat
    tuple val(id), path("*.phased.vcf.gz*"), emit: phasedvcf
    tuple val(id), path("*.phase.report"), emit: report 

    tag "$id"
    publishDir "${params.outdir}/$id/phase/", mode: 'link'

    // publishDir (
    //     path: "${params.outdir}/$id/phase/", 
    //     saveAs: { fn ->
    //         if (fn.contains("vcf.gz") || fn.contains("hapblock") || fn.contains("hapcut_stat")) {"$fn"}
    //         else {"../../report/$id/$fn"}
    //     }
    // )
    // cache false
    script:
    def prefix = "${id}.${params.align_tool}.${params.var_tool}"
    def fai = params.ref.startsWith('/') ? "${params.ref}.fai" : "${params.DB}/${params.ref}/reference/${params.ref}.fa.fai"
    def chr1 = (params.ref == "hs37d5") ? "" : "chr"
    def py = "${params.SCRIPT}/calculate_haplotype_statistics_CWX.py -h1 \$hapblocks -v1 \$pvcfs -v2 \$pvs --indels >> ${prefix}.hapcut_stat.txt"
    """
    lfs=""
    hapblocks=""
    stat2s=""
    pvs=""
    pvcfs=""
    vcfs=""

    if [ "${params.chr}" = "all" ]; then
        for i in {1..22} X;do
            lfs="\$lfs ${prefix}.${chr1}\${i}.lf"
            hapblocks="\$hapblocks ${prefix}.${chr1}\${i}.hapblock"
            stat2s="\$stat2s ${prefix}.${chr1}\${i}.hapcut_stat.txt"
            pvs="\$pvs ${params.DB}/hg38/phasedvcf/hg38.${params.std}.${chr1}\${i}.vcf.gz"
            vcfs="\$vcfs ${prefix}.${chr1}\${i}.vcf.gz"
            pvcfs="\$pvcfs ${prefix}.${chr1}\${i}.hapblock.phased.VCF.gz"
        done
    else 
        lfs="${prefix}.${params.chr}.lf"
        hapblocks="${prefix}.${params.chr}.hapblock"
        stat2s="${prefix}.${params.chr}.hapcut_stat.txt"
        pvs="${params.DB}/${params.ref}/phasedvcf/${params.ref}.${params.std}.${params.chr}.vcf.gz"
        vcfs="${prefix}.${params.chr}.vcf.gz"
        pvcfs="${prefix}.${params.chr}.hapblock.phased.VCF.gz"
    fi

    cat \$lfs > ${prefix}.lf
    cat \$hapblocks > ${prefix}.hapblock
    cat \$stat2s > ${prefix}.hapcut_stat.txt

    echo "combine all chrs" >> ${prefix}.hapcut_stat.txt

    python3 $py

 	${params.BIN}bcftools concat *phased.VCF.gz -O b -o tmp.vcf.gz 
    zcat tmp.vcf.gz | grep '^#' > header
    zcat tmp.vcf.gz | grep -v '^#' | sort -k1,1d -k2,2n > body
    cat header body |bgzip -c > ${prefix}.phased.vcf.gz
    rm tmp.vcf.gz
    
    ${params.BIN}tabix ${prefix}.phased.vcf.gz

    
    python3 ${params.SCRIPT}/phase.py ${prefix} $fai > ${prefix}.phase.report
    """
}
process phaseCatRef {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    
    input:
    path(txt)
    tuple val(id), path(vcfs), path(pvs), path(lfs), path(hbs)

    output:
    tuple val(id), path("*hapblock"), emit: hb
    tuple val(id), path("*.phased.vcf.gz*"), emit: phasedvcf
    tuple val(id), path("*.phase.report"), emit: report 

    tag "$id"
    publishDir "${params.outdir}/$id/phase/", mode: 'link'

    script:
    def prefix = "${id}.${aligner}.${varcaller}"
    """
    lfs=""
    hapblocks=""
    pvs=""
    pvcfs=""
    vcfs=""

    while IFS= read -r i; do
        lfs="\$lfs ${prefix}.\${i}.lf"
        hapblocks="\$hapblocks ${prefix}.\${i}.hapblock"
        vcfs="\$vcfs ${prefix}.\${i}.vcf.gz"
        pvcfs="\$pvcfs ${prefix}.\${i}.hapblock.phased.VCF.gz"
    done < $txt


    cat \$lfs > ${prefix}.lf
    cat \$hapblocks > ${prefix}.hapblock


 	${params.BIN}bcftools concat *phased.VCF.gz -O b -o tmp.vcf.gz 
    zcat tmp.vcf.gz | grep '^#' > header
    zcat tmp.vcf.gz | grep -v '^#' | sort -k1,1d -k2,2n > body
    cat header body |bgzip -c > ${prefix}.phased.vcf.gz
    rm tmp.vcf.gz
    
    ${params.BIN}tabix ${prefix}.phased.vcf.gz

    python3 ${params.SCRIPT}/calc_n50.py -h1 \$hapblocks -v1 \$pvcfs --indels >> ${prefix}.hapcut_stat.txt
    python3 ${params.SCRIPT}/phase.py ${prefix} 0 > ${prefix}.phase.report # allhetsnp, phasedhetsnp, allhetindel, phasedhetindel, fbn50, fbnum
    """
}
