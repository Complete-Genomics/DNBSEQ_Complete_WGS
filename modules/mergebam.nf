workflow WF_mergebam {
    take:
    ch_bams

    main:
    if (params.just_combine) {
        combinebam(ch_bams).set {ch_mergebam}
    } else {
        intersect(ch_bams).set {ch_bed}
        mergeBam(ch_bams.join(ch_bed)).set {ch_mergebam}
    }
    

    emit:
    ch_mergebam
}

process intersect {
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(stlfrbam), path(pfbam)

    output:
    tuple val(id), path("*.intersect.bed") //demo.stlfr.bwa.cov10.bed

    tag "$id"
    publishDir "${params.outdir}/$id/align/", mode: 'link'
    
    script:
    def aligner = params.align_tool
    stlfrbam = stlfrbam.first() //demo.stlfr.bwa.bam
    pfbam = pfbam.first()
    """
    ${params.BIN}samtools depth -@ ${task.cpus} $stlfrbam $pfbam | awk -v cov="${params.PF_lt_stLFR_depth}" '\$3 >= cov && \$4 < cov {print \$1"\\t"\$2"\\t"\$2}' | sort -k1,1 -k2,2n | \\
    ${params.BIN}bedtools merge > ${id}.${aligner}.cov${params.PF_lt_stLFR_depth}.intersect.bed

    """
    stub:
    "touch ${id}.${params.align_tool}.cov${params.PF_lt_stLFR_depth}.intersect.bed"
}
process mergeBam {
    cpus params.cpu2
    memory params.MEM3 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(stlfrbam), path(pfbam), path(bed) 

    output:
    tuple val(id), path("*.merge.bam*") 
    tag "$id"

    publishDir "${params.outdir}/$id/align/", mode: 'link'

    script:
    def aligner = params.align_tool
    def pfbam = pfbam.first()
    def stlfrbam = stlfrbam.first()
    """
    ${params.BIN}samtools view -hb -L $bed $stlfrbam > lfr_lfr10_pf10.bam && \\
    ${params.BIN}samtools index lfr_lfr10_pf10.bam && \\
    ${params.BIN}samtools addreplacerg -w -O BAM -r '@RG\\tID:${id}\\tSM:sample' -o new.bam lfr_lfr10_pf10.bam && \\
    ${params.BIN}samtools index new.bam && \\
    ${params.BIN}samtools reheader $pfbam new.bam > new2.bam && \\
    ${params.BIN}samtools index new2.bam && \\
    ${params.BIN}samtools addreplacerg -O BAM -r '@RG\\tID:sample\\tSM:sample' -o new3.bam new2.bam && \\
    ${params.BIN}samtools index new3.bam && \\
    ${params.BIN}samtools merge -@ ${task.cpus} -f ${id}.${aligner}.merge.bam $pfbam new3.bam && \\
    ${params.BIN}samtools index -@ ${task.cpus} ${id}.${aligner}.merge.bam

    """
    stub:
    "touch ${id}.${params.align_tool}.merge.bam ${id}.${params.align_tool}.merge.bam.bai"
}

process combinebam {
    cpus params.cpu2
    memory params.MEM3 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(stlfrbam), path(pfbam)

    output:
    tuple val(id), path("${id}.lariat.merge.bam*") 

    publishDir "${params.outdir}/$id/align/", mode: 'link'

    script:
    def pfbam = pfbam.first()
    def stlfrbam = stlfrbam.first()
    def aligner = params.align_tool
    """
    ${params.BIN}samtools merge -@ ${task.cpus} -f ${id}.${aligner}.merge.bam $pfbam $stlfrbam && \\
    ${params.BIN}samtools index -@ ${task.cpus} ${id}.${aligner}.merge.bam

    """
    stub:
    "touch ${id}.lariat.merge.bam ${id}.lariat.merge.bam.bai"
}
