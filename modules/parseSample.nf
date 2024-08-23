import java.nio.file.Paths

workflow parse_sample {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    toCsv(samplesheet).set {ch_x}
    // Channel.fromPath( "${params.outdir}/samplesheet.csv" )
    // Channel.fromPath(samplesheet)
        // .view()
    ch_x
        .splitCsv ( header:true, sep:',' ) // dict: [sample:hg002, stlfr1: path, ...]
        .map { create_fastq_channel(it) }   //       [test1, [stlfr1, stlfr2, pf1, pf2]]
        .set { reads }                      
    emit: reads   
}
workflow parse_sample_frombam {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    toCsv(samplesheet).set {ch_x}

    // Channel.fromPath( "${params.outdir}/samplesheet.csv" )
    // Channel.fromPath(samplesheet)
        // .view()
    ch_x
        .splitCsv ( header:true, sep:',' ) // dict: [sample:hg002, stlfr1: path, ...]
        .map { create_channel_frombam(it) }   //       [test1, [stlfr1, stlfr2, pf1, pf2]]
        .set { reads }                      
    emit: reads   
}

process tosamplelist {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

    input:
    path(samplesheet)

    output:
    path("input.samplesheet")
    // publishDir "${params.outdir}", mode: "copy"
    // cache false
    "${params.BIN}python ${params.SCRIPT}/tosamplelist.py $samplesheet input.samplesheet"
}
process toCsv {
    executor = 'local'
    container false
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

    input:
    path(samplesheet)

    output:
    path("samplesheet.csv")
    // publishDir "${params.outdir}", mode: "copy"

    "sed 's/[[:blank:]]\\+/,/g' $samplesheet > samplesheet.csv"
}
process bam {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

    input:
    tuple val(id), path(stlfrbam), path(pfbam)

    output:
	tuple val(id), path("${id}_stlfr.bam*"),emit: stlfr
	tuple val(id), path("${id}_pf.bam*"),emit: pf

	tag "$id"

    // publishDir "${params.outdir}", mode: "copy"

    script:
	"""
	mv $stlfrbam ${id}_stlfr.bam
	mv $pfbam ${id}_pf.bam

    ${params.BIN}samtools index ${id}_stlfr.bam
    ${params.BIN}samtools index ${id}_pf.bam
	"""
}
def create_fastq_channel(LinkedHashMap row) {

    def resolvePath = { path ->
        def p = Paths.get(path)
        return p.isAbsolute() ? p : workflow.launchDir.parent.resolve(p).toAbsolutePath()
    }

    def stlfr1 = resolvePath(row.stlfr1).toString()
    def stlfr2 = resolvePath(row.stlfr2).toString()
    def pcrfree1 = resolvePath(row.pcrfree1).toString()
    def pcrfree2 = resolvePath(row.pcrfree2).toString()

    return [row.sample, stlfr1, stlfr2, pcrfree1, pcrfree2]
}

def create_channel_frombam(LinkedHashMap row) {
    def resolvePath = { path ->
        def p = Paths.get(path)
        return p.isAbsolute() ? p : workflow.launchDir.parent.resolve(p).toAbsolutePath()
    }
    def stlfrbam = resolvePath(row.stlfrbam).toString()
    def pfbam = resolvePath(row.pfbam).toString()
    return [row.sample, stlfrbam, pfbam]
}