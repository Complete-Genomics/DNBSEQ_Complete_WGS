process bam2bed {
    
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path("bam.bed")

    script:
    def bam = bam.first()
    """
    ${params.BIN}bedtools bamtobed -i $bam > bam.bed
    """
}

process coverage {
    
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(lib)
    tuple val(id), path(bam)

    output:
    tuple val(id), path("${id}.*.bed")

    tag "$id, $lib"
    publishDir "${params.outdir}/$id/align/", mode:'link'

    script:
    bam = bam.first()
    def fai = "${params.DB}/${params.ref}/reference/${params.ref}.fa.fai"
    def bed = "${params.DB}/hg38/GRCh38_CMRG_benchmark_gene_coordinates.bed"
    """
    ${params.BIN}bedtools coverage -sorted -g $fai -a $bed -b $bam -hist > ${id}.${lib}.cmrg.hist.bed 
    """
    stub:
    "touch ${id}.${lib}.cmrg.hist.bed "
}

process coverageMean {
    
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(lib)
    tuple val(id), path(bam)

    output:
    tuple val(id), path("${id}.*.bed")

    tag "$id"
    publishDir "${params.outdir}/$id/align/", mode:'link'

    script:
    bam = bam.first()
    def fai = "${params.DB}/${params.ref}/reference/${params.ref}.fa.fai"
    def bed = "${params.DB}/hg38/GRCh38_CMRG_benchmark_gene_coordinates.bed"
    """
    ${params.BIN}bedtools coverage -sorted -g $fai -a $bed -b $bam -mean > ${id}.${lib}.cmrg.mean.bed 
    """
    stub:
    "touch ${id}.${lib}.cmrg.mean.bed "
}

process coverageAvg {
    
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    tuple val(id), path(pfhistbed), path(pfmeanbed), path(mergehistbed), path(mergemeanbed)

    output:
    tuple val(id), path("${id}.cmrg.cov")

    tag "$id"
    publishDir "${params.outdir}/$id/align/", mode:'link'

    script:
    """
    ${params.BIN}python3 ${params.SCRIPT}/cmrg.py $pfhistbed $mergehistbed $id 

    pfavgcov=`awk '{sum+=\$2;n++}END{print sum/n}' ${id}.cmrg.pf.cov`
    mergeavgcov=`awk '{sum+=\$2;n++}END{print sum/n}' ${id}.cmrg.merge.cov`

    pfavgdepth=`awk '{sum += \$NF; count++} END {if (count > 0) print sum / count}' $pfmeanbed`
    mergeavgdepth=`awk '{sum += \$NF; count++} END {if (count > 0) print sum / count}' $mergemeanbed`

    echo -e "\$pfavgcov\\t\$mergeavgcov\\t\$pfavgdepth\\t\$mergeavgdepth" > ${id}.cmrg.cov
        
    #${params.BIN}samtools depth -b cov0.7.bed bam  > depth
    #python3 ${params.SCRIPT}/genedepth.py \$genebed depth > out
    """
    stub:
    "touch ${id}.cmrg.cov "
}
process mosdepth {
    cpus params.cpu2
    memory params.MEM2 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(lib)
    tuple val(id), path(bam)

    output:
    tuple val(id), path("*.regions.bed.gz")

    tag "$id"

    publishDir "${params.outdir}/$id/align/", mode:'link'

    script:
    bam = bam.first()
    def prefix = "${id}.${lib}"
    def bed = "${params.DB}/${params.ref}/merge_1216_${params.ref}.bed"
    """
    /hwfssz8/MGI_BIOINFO/USER/biancai/pwd/envs/mosdepth/bin/mosdepth -t ${task.cpus} -n --by $bed $prefix $bam 2> ${prefix}.mosdepth.log
    ${params.BIN}Rscript ${params.SCRIPT}/1.visualization_density_plot_coverage.R ${prefix}.regions.bed.gz ${prefix}.mosdepth.pdf
    """
}
process mosdepthPlot {
    cpus params.cpu2
    memory params.MEM2 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
    input:
    val(lib)
    tuple val(id), path(bed1), path(bed2)

    output:
    tuple val(id), path("${id}.${lib}.mosdepth.pdf")

    tag "$id"

    publishDir "${params.outdir}/report/$id/", mode:'link'

    script:
    bam = bam.first()
    def prefix = "${id}.${lib}"
    def bed = "${params.DB}/${params.ref}/merge_1216_${params.ref}.bed"
    """
    ls *regions.bed.gz > paths
    ${params.BIN}Rscript ${params.SCRIPT}/1.visualization_density_plot_coverage.R paths ${id}.${lib}.mosdepth.pdf
    """
}