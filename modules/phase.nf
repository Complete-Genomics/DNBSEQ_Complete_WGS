process getchrs {
    executor = 'local'
	container false

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
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    

    input:
    val(aligner)
    val(varcaller)
    tuple val(id), val(chr), path(bam), path(vcf) //demo.stlfr.gatk.chr1.vcf.gz

    output:
    tuple val(id), path("*.VCF.gz"), emit: phasedvcf
    tuple val(id), path("*.lf"), emit: lf
    tuple val(id), path("*.hapblock"), emit: hapblock
    tuple val(id), val(chr), path(bam), path("*.hapblock"), emit: svpre
    tuple val(id), path("*.hapcut_stat.txt"), emit: stat

    tag "$id, $aligner, $varcaller, $chr"
    // publishDir "${params.outdir}/$id/phase/phasesplit/"

    // publishDir (
    //     path: "${params.outdir}/$id/phase/", 
    //     saveAs: { fn ->
    //         if (fn.contains("lf") || fn.contains("hapblock") || fn.contains("hapcut_stat")) {"phasesplit/$fn"}
    //         else {"svsplit/$fn"}
    //     }
    // )
    
    // cache false

    script:
    def bam = bam.first()
    def prefix = "${id}.${aligner}.${varcaller}.${chr}"
    def pv = "${params.DB}/${params.ref}/phasedvcf/${params.ref}.${params.std}.${chr}.vcf.gz"
    def fa = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    def fai = "${fa}.fai"
    cmd = """
    #export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:${params.DB}/htslib
    ${params.BIN}gzip -dc $vcf > tmp.vcf
    ${params.BIN_HAPCUT2}extractHAIRS --indels 1 --10X 1 --bam $bam --VCF tmp.vcf --out ${prefix}.unlinked_frag # --triallelic 1

    python3 ${params.SCRIPT_HapCUT2}LinkFragments.py \\
        --bam $bam --VCF tmp.vcf --fragments ${prefix}.unlinked_frag --out ${prefix}.lf -d 100000

    ${params.BIN_HAPCUT2}HAPCUT2 --nf 1 --fragments ${prefix}.lf --vcf tmp.vcf --output ${prefix}.hapblock

    ${params.BIN}bgzip ${prefix}.hapblock.phased.VCF
    rm tmp.vcf
    echo $chr > ${prefix}.hapcut_stat.txt
    """
    if (params.stLFR_only) {
        cmd += """
        python3 ${params.SCRIPT}/calculate_haplotype_statistics.reseq.py \\
        -h1 ${prefix}.hapblock -v1 ${prefix}.hapblock.phased.VCF.gz -f1 ${prefix}.lf -pv $pv -c $fai >> ${prefix}.hapcut_stat.txt
        """
    } else {
        cmd += """
        if [ -s ${prefix}.hapblock ]; then
            python3 ${params.SCRIPT}/calculate_haplotype_statistics_CWX.py \\
            -h1 ${prefix}.hapblock -v1 ${prefix}.hapblock.phased.VCF.gz -v2 $pv --indels >> ${prefix}.hapcut_stat.txt
        fi
        """
    }
    return cmd
    stub:
    "touch ${id}.${aligner}.${varcaller}.${chr}.VCF.gz ${id}.${aligner}.${varcaller}.${chr}.lf ${id}.${aligner}.${varcaller}.${chr}.hapblock ${id}.${aligner}.${varcaller}.${chr}.hapcut_stat.txt"
}

process eachstat_phase {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    

    input:
    // val(aligner)
    // val(varcaller)
    tuple val(id), path(hapcutstat), path(vcfs)

    output:
    tuple val(id), path("*.haplotype.xls")

    // tag "$id, $aligner, $varcaller"
    tag "$id"
    publishDir "${params.outdir}/report/$id/", mode: 'link'

    // cache false

    script:
    // def prefix = "${id}.${aligner}.${varcaller}"
    """    
    # summarize phasing data from id.hapcut_stat.txt
    ${params.BIN}python ${params.SCRIPT}/stat/eachstat_phase.py $id $hapcutstat > ${id}.haplotype.xls
    """
}
process hapKaryotype {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    

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
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
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
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    

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
process phase_cat {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(aligner)
    val(varcaller)
    tuple val(id), path(vcfs), path(pvs), path(lfs), path(hbs), path(stats)

    output:
    tuple val(id), path("*hapblock"), emit: hb
    tuple val(id), path("${id}.${aligner}.${varcaller}.hapcut_stat.txt"), emit: hapcutstat
    path("*.phased.vcf.gz*")
    tuple val(id), path("*.phase.report"), emit: report 

    tag "$id, $aligner, $varcaller"
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
    def prefix = "${id}.${aligner}.${varcaller}"
    def fai = "${params.DB}/${params.ref}/reference/${params.ref}.fa.fai"
    def fa = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    def chr1 = (params.ref == "hs37d5") ? "" : "chr"
    def py = (params.stLFR_only) ? "${params.SCRIPT}/calculate_haplotype_statistics.reseq.py -h1 \$hapblocks -v1 \$pvcfs -f1 \$lfs -pv \$pvs -c $fai >> ${prefix}.hapcut_stat.txt" : "${params.SCRIPT}/calculate_haplotype_statistics_CWX.py -h1 \$hapblocks -v1 \$pvcfs -v2 \$pvs --indels >> ${prefix}.hapcut_stat.txt"
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
            pvs="\$pvs ${params.DB}/${params.ref}/phasedvcf/${params.ref}.${params.std}.${chr1}\${i}.vcf.gz"
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
    stub:
    "touch *.phase.report *.phased.vcf.gz ${id}.${aligner}.${varcaller}.VCF.gz ${id}.${aligner}.${varcaller}.lf ${id}.${aligner}.${varcaller}.hapblock ${id}.${aligner}.${varcaller}.hapcut_stat.txt"
}
process phaseCatRef {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(aligner)
    val(varcaller)
    path(txt)
    tuple val(id), path(vcfs), path(pvs), path(lfs), path(hbs)

    output:
    tuple val(id), path("*hapblock"), emit: hb
    path("*.phased.vcf.gz*")
    // tuple val(id), path("*.phase.report"), emit: report 

    tag "$id, $aligner, $varcaller"
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

    # python3 ${params.SCRIPT}/phase.py ${prefix} \$fai > ${prefix}.phase.report
    """
    stub:
    "touch *.phase.report *.phased.vcf.gz ${id}.${aligner}.${varcaller}.VCF.gz ${id}.${aligner}.${varcaller}.lf ${id}.${aligner}.${varcaller}.hapblock ${id}.${aligner}.${varcaller}.hapcut_stat.txt"
}
process phaseCat_cwx {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(aligner)
    val(varcaller)
    tuple val(id), path(vcfs), path(pvs), path(lfs), path(hbs), path(stats)

    // output:
    // tuple val(id), path("*hapblock"), emit: hb
    // tuple val(id), path("${id}.${aligner}.${varcaller}.hapcut_stat.txt"), emit: hapcutstat
    // path("*.phased.vcf.gz*")
    // tuple val(id), path("*.phase.report"), emit: report 

    // tag "$id, $aligner, $varcaller"
    // publishDir "${params.outdir}/$id/phase/", mode: 'link'

    // publishDir (
    //     path: "${params.outdir}/$id/phase/", 
    //     saveAs: { fn ->
    //         if (fn.contains("vcf.gz") || fn.contains("hapblock") || fn.contains("hapcut_stat")) {"$fn"}
    //         else {"../../report/$id/$fn"}
    //     }
    // )

    script:
    def prefix = "${id}.${aligner}.${varcaller}"
    def fai = "${params.DB}/${params.ref}/reference/${params.ref}.fa.fai"
    def fa = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    def chr1 = (params.ref == "hs37d5") ? "" : "chr"
    """
    lfs=""
    hapblocks=""
    stat2s=""
    pvs=""
    pvcfs=""
    vcfs=""

    for i in {1..22} X;do
        lfs="\$lfs ${prefix}.${chr1}\${i}.lf"
        hapblocks="\$hapblocks ${prefix}.${chr1}\${i}.hapblock"
        stat2s="\$stat2s ${prefix}.${chr1}\${i}.hapcut_stat.txt"
        pvs="\$pvs ${params.DB}/${params.ref}/phasedvcf/${params.ref}.${params.std}.${chr1}\${i}.vcf.gz"
        vcfs="\$vcfs ${prefix}.${chr1}\${i}.vcf"
        pvcfs="\$vcfs ${prefix}.${chr1}\${i}.hapblock.phased.VCF"
    done

    cat \$lfs > ${prefix}.lf
    cat \$hapblocks > ${prefix}.hapblock
    cat \$stat2s > ${prefix}.hapcut_stat.txt

    echo "combine all chrs" >> ${prefix}.hapcut_stat.txt

    python3 ${params.SCRIPT}/calculate_haplotype_statistics_CWX.py \\
        -h1 \$hapblocks -v1 \$pvcfs -v2 \$pvs --indels >> ${prefix}.hapcut_stat.txt

    """
}