import java.nio.file.Paths

workflow parse_sample {
    take:
    samplesheet

    main:
    toCsv(samplesheet)
        .splitCsv(header:true)
        .flatMap { row ->               // flatMap：一行变多行
            def id = row.sample
            def out = []

            // stlfr2 column is ambiguous:
            //   - row has stlfr1 + stlfr2 → regular PE stLFR (stlfr2 = R2)
            //   - row has stlfr2 only (no stlfr1) → SE stLFR2 library
            def hasStlfr1 = row.containsKey('stlfr1') && row['stlfr1']?.trim()

            // 先按列名把文件分组
            def byLibType = [:]         // [stlfr:[fq:[], bam:[]], stlfr2:[fq:[], bam:[]], pf:[fq:[], bam:[]]]
            row.each { hdr, pathStr ->
                if (!pathStr) return     // 空路径跳过
                def trimmed = pathStr.trim()
                if (!trimmed) return
                def lib
                if ((hdr == 'stlfr2' && !hasStlfr1) || hdr.startsWith('stlfr2') && hdr != 'stlfr2') {
                    lib = 'stlfr2'   // SE stLFR2 (stlfr2 alone) or PE SE600 (stlfr21/stlfr22)
                } else if (hdr.contains('stlfr')) {
                    lib = 'stlfr'    // regular PE stLFR (stlfr1=R1, stlfr2=R2)
                } else if (hdr.startsWith('pf') || hdr.startsWith('pcrfree')) {
                    lib = 'pf'
                } else {
                    lib = null
                }
                if (!lib) return
                def type = hdr.endsWith('bam') ? 'bam' : 'fq'

                if (!byLibType[lib])      byLibType[lib] = [:]
                if (!byLibType[lib][type]) byLibType[lib][type] = []

                byLibType[lib][type] << [hdr: hdr, path: file(trimmed)]
            }

            // 遍历 (lib,type) 生成独立记录
            byLibType.each { lib, typeMap ->
                typeMap.each { type, files ->
                    def meta = [id: id, lib: lib, type: type]

                    if (type == 'fq') {
                        // fastq：按列名最后一位 1/2 排序，保证 R1/R2 顺序
                        def ordered = files.sort { it.hdr[-1] as int }.collect{ it.path }
                        out << [meta, ordered]          // 列表
                    } else {
                        def bamFile = files[0].path
                        def bamStr  = bamFile.toString()
                        def baiFile = file("${bamStr}.bai")

                        out << [meta, [bamFile, baiFile]]
                    }
                }
            }
            return out
        }
        .set { data }

    emit:
    data
}

// workflow parse_sample {
//     take:
//     samplesheet // file: /path/to/samplesheet.csv

//     main:
//     toCsv(samplesheet)
//         .splitCsv ( header:true ) // dict: [sample:hg002, stlfr1: path, ...]
//         .map { func(it) }   //       
//         .set { reads }  

    
//     emit: 
//     reads  
// }

workflow parse_sample_frombam {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    toCsv(samplesheet).set {ch_x}

    ch_x
        .splitCsv(header:true, sep:',')
        .multiMap { row ->
            bam:    create_channel_frombam(row)
            fq_pf:  create_channel_frombam_fq(row)   // pf1/pf2 columns, may be null
        }
        .set { parsed }

    parsed.bam.set { bam }
    parsed.fq_pf
        .filter { id, fqs -> fqs != null }
        .set { fq_pf }

    emit:
    bam   = bam
    fq_pf = fq_pf
}

workflow bam {
    take:
    ch_bam

    main:
    ch_bam.map { id, stlfrbam, pfbam ->
        [id, [file(stlfrbam), file("${stlfrbam}.bai")]]
    }.set { stlfr }

    ch_bam.map { id, stlfrbam, pfbam ->
        [id, [file(pfbam), file("${pfbam}.bai")]]
    }.set { pf }

    emit:
    stlfr = stlfr
    pf    = pf
}

workflow bam2 {
    take:
    ch_bam

    main:
    ch_bam.map { id, stlfrbam, mergebam ->
        [id, [file(stlfrbam), file("${stlfrbam}.bai")]]
    }.set { stlfr }

    ch_bam.map { id, stlfrbam, mergebam ->
        [id, [file(mergebam), file("${mergebam}.bai")]]
    }.set { merge }

    emit:
    stlfr = stlfr
    merge = merge
}

process tosamplelist {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    path(samplesheet)

    output:
    path("input.samplesheet")
    // publishDir "${params.outdir}", mode: "copy"
    // cache false
    "${params.BIN}python ${params.SCRIPT}/tosamplelist.py $samplesheet input.samplesheet"
}
process toCsv {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    path(samplesheet)

    output:
    path("samplesheet.csv")
    // publishDir "${params.outdir}", mode: "copy"

    "sed 's/[[:blank:]]\\+/,/g' $samplesheet > samplesheet.csv"
}

def func(LinkedHashMap row) {
    def resolvePath = { path ->
        def p = Paths.get(path)
        return p.isAbsolute() ? p : workflow.launchDir.parent.resolve(p).toAbsolutePath()
    }
    def newRow = [:]
    row.each { key, value ->
        if (key == 'sample') {
            newRow[key] = value
        } else {
            // 处理相对路径（如 ../demo.fq）
            path = resolvePath(value).toString()

            if (!file(path).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> ${path} does not exist!\n"
            }

            newRow[key] = path
        }
    }
    return newRow
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

    if (!file(stlfr1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> ${stlfr1} does not exist!\n"
    }
    if (!file(stlfr2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> ${stlfr2} does not exist!\n"
    }
    if (!file(pcrfree1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> ${pcrfree1} does not exist!\n"
    }
    if (!file(pcrfree2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> ${pcrfree2} does not exist!\n"
    }
    return [row.sample, stlfr1, stlfr2, pcrfree1, pcrfree2]
}

def create_channel_frombam_fq(LinkedHashMap row) {
    def resolvePath = { path ->
        def p = java.nio.file.Paths.get(path)
        return p.isAbsolute() ? p : workflow.launchDir.parent.resolve(p).toAbsolutePath()
    }
    def pf1str = row['pcrfree1']?.trim()
    def pf2str = row['pcrfree2']?.trim()
    if (!pf1str || !pf2str) return [row.sample, null]

    def pf1 = file(resolvePath(pf1str).toString())
    def pf2 = file(resolvePath(pf2str).toString())
    if (!pf1.exists()) exit 1, "ERROR: pcrfree1 not found: ${pf1}\n"
    if (!pf2.exists()) exit 1, "ERROR: pcrfree2 not found: ${pf2}\n"
    return [row.sample, [pf1, pf2]]
}

def create_channel_frombam(LinkedHashMap row) {
    def resolvePath = { path ->
        def p = Paths.get(path)
        return p.isAbsolute() ? p : workflow.launchDir.parent.resolve(p).toAbsolutePath()
    }
    if (!params.fromMergedBam) {
        def stlfrbam = resolvePath(row.stlfrbam).toString()
        def pfbam = resolvePath(row.pfbam).toString()

        if (!file(stlfrbam).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> ${stlfrbam} does not exist!\n"
        }
        if (!file(stlfrbam).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> ${stlfrbam} does not exist!\n"
        }

        return [row.sample, stlfrbam, pfbam]

    } else { // fromMergedBam
        def stlfrbam = resolvePath(row.stlfrbam).toString()
        def mergebam = resolvePath(row.mergebam).toString()
        if (!file(stlfrbam).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> ${stlfrbam} does not exist!\n"
        }
        if (!file(mergebam).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> ${mergebam} does not exist!\n"
        }
        return [row.sample, stlfrbam, mergebam]
    }
}
