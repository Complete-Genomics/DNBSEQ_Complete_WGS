include {
    qc as qc_pf;
    qc_stlfr_stats                                 } from "${params.MOD}/qc"
include {
    readLen as readLenPf;
    basecount as basecountPf;
    splitfq;
    samplePfFq;
    sampleStlfrFq       } from "${params.MOD}/fq"
include { mapq          } from "${params.MOD}/bam"


workflow WF_align_stlfr {
    take:
    ch_stlfrfq

    main:
    // ch_stlfrfq.map { meta, path ->
    //     return [meta.id, path]
    // }.set {ch_stlfrfq}
    needSplit(ch_stlfrfq).branch {id, reads, flag ->
        split: flag == '1'
            [id, reads, true]
        no_split: flag == '0'
            [id, reads, false]
    }.set {ch_dosplit}

    barcode_split(ch_dosplit.split)
    splitRate(ch_dosplit.no_split)

    ch_splitfq  = barcode_split.out.reads   .mix(splitRate.out.reads)
    splitLog    = barcode_split.out.log     .mix(splitRate.out.log)

    if (params.sampleFq) {
        qc_stlfr_stats(ch_splitfq).basecnt.set {ch_stlfrbasecount}
        qc_stlfr_stats.out.rlen.set {ch_stLFRreadLen}
        sampleStlfrFq(ch_stlfrbasecount.join(ch_stLFRreadLen).join(ch_splitfq)).set {ch_splitfq} 
    }
    if (params.align_tool == 'lariat') {
        if (params.lariatStLFRBC) { //use stLFR bc for lariat
            splitfq(ch_splitfq).fq1s.transpose().set {ch_fq1s} //[id, num, fq1, fq2]
            splitfq.out.fq2s.transpose().set {ch_fq2s}
            ch_fq1s.join(ch_fq2s).set {ch_splitstlfrfq}
            lariatBC(ch_splitstlfrfq).groupTuple().set {ch_splitlariatfqs}
            mergeFq(ch_splitlariatfqs).set {ch_lariatfq}
            // ch_splitlariatfqs.map {it -> it[1]}.collect().view()
        } else {
            tofake10xHash(splitLog).set {ch_hash}
            if (params.lariatSplitFqNum != 1) { 
                splitfq(ch_splitfq).fq1s.transpose().set {ch_fq1s} //[id, num, fq1, fq2]
                splitfq.out.fq2s.transpose().set {ch_fq2s}
                ch_fq1s.join(ch_fq2s).set {ch_splitstlfrfq}
                // ch_splitstlfrfq.view()
                tofake10x(ch_splitstlfrfq.combine(ch_hash, by:0)).reads.set {ch_splitfake10xfq}
                fake10x2lariat(ch_splitfake10xfq).groupTuple().set {ch_splitlariatfqs}
                mergeFq(ch_splitlariatfqs).set {ch_lariatfq}
            } else { fake10x2lariat(tofake10x(ch_splitfq.join(ch_hash))).set {ch_lariatfq} }
        }
        sortbam(lariat(ch_lariatfq)).set {ch_lariatbam0} 
        markdup('stlfr', 'lariat', ch_lariatbam0).set {ch_stlfrbam}
    } else { // bwa
        if (params.use_megabolt) {
            bwaMegabolt('stlfr', ch_splitfq).set {ch_stlfrBwaBam}
        } else {
            bwa('stlfr', ch_splitfq).set {ch_stlfrBwaBam}
            markdup('stlfr', 'bwa', ch_stlfrBwaBam).set {ch_stlfrbam}
        }
    }
    
    if (params.sampleBam) { sampleBamStlfrLariat('stlfr', 'lariat', ch_stlfrbam).set {ch_stlfrbam} } 

    emit:
    ch_stlfrbam
}
workflow WF_align_pf {
    take:
    ch_pffq

    main:                                            
    // ch_pffq.map { meta, path ->
    //     return [meta.id, path]
    // }.set {ch_pffq}

    qc_pf('pf', ch_pffq).reads.set {ch_qcpffq} 

    if (params.sampleFq) { 
        qc_pf.out.bssq.set {ch_pfbssq}
        readLenPf(ch_pfbssq).set {ch_PFreadLen} 
        basecountPf(ch_pfbssq).set {ch_pfbasecount}
        samplePfFq(ch_pfbasecount.join(ch_PFreadLen).join(ch_qcpffq)).set {ch_pffq}
    }

    if (params.pfAligner == 'bwa') {
        if (params.use_megabolt) {
            bwaMegaboltPf('pf', ch_pffq).set {ch_pfbam}
        } else {
            bwa('pf', ch_pffq).set {ch_pfsortbam}
            markdup('pf', 'bwa', ch_pfsortbam).set {ch_pfbam} 
        } 
    } else if (params.pfAligner == 'vg') {
        if (params.ref != 'hg38' && !params.ref.contains('GRCh38')) { exit 1, 'graph aligner only support hg38 ref!'}
        ch_pffq.map { id, reads -> [id, reads, true] }.set {ch_pffq_typed}
        kff(ch_pffq_typed.map { id, reads, is_pe -> [id, reads] }).set {ch_kff}
        vg(ch_kff.join(ch_pffq_typed)).set {ch_pfbam}
        markdup('pf', 'vg', ch_pfbam).set {ch_pfbam}
    } else {
        exit 1, 'pf aligner only supports bwa and vg!'
    }
    if (params.sampleBam) { sampleBam('pf', 'bwa', ch_pfbam).set {ch_pfbam} }
    mapq(ch_pfbam).set {ch_pfbam}

    emit:
    ch_pfbam
}
// stLFR single-end 600/700bp reads
// - fastq columns in samp.list: stlfr21 / stlfr22  (same PE format as regular stLFR)
// - barcodes are at the 3' end of each read
// - Lariat cannot handle SE → align with vg giraffe (same as PF)
// - vg process already strips GRCh38#0# from BAM header (sed 's/GRCh38#0#//g')
workflow WF_align_stlfr2 {
    take:
    ch_stlfr2fq   // [id, [r1]] — single-end raw FASTQ (barcodes embedded in sequence)

    main:
    if (params.ref != 'hg38' && !params.ref.contains('GRCh38')) {
        exit 1, 'stlfr2 aligner only supports hg38/GRCh38 ref (vg giraffe)!'
    }

    // SE600: extract barcodes from read sequence, encode in read name (#bc1_bc2_bc3)
    barcode_split_se600(ch_stlfr2fq).set { ch_stlfr2_split }

    // Pass split FASTQ as single-element list to avoid Nextflow path-staging collision
    ch_stlfr2_split.map { id, fq -> [id, [fq], false] }
                   .set { ch_stlfr2fq_typed }   // [id, [fq], is_pe=false]

    kff(ch_stlfr2fq_typed.map { id, fqs, is_pe -> [id, fqs] }).set { ch_kff2 }
    vg(ch_kff2.join(ch_stlfr2fq_typed)).set { ch_stlfr2bam0 }
    addBxSe600(ch_stlfr2bam0).set { ch_stlfr2bam0 }
    markdup('stlfr2', 'vg', ch_stlfr2bam0).set { ch_stlfr2bam }

    if (params.sampleBam) { sampleBam('stlfr2', 'vg', ch_stlfr2bam).set { ch_stlfr2bam } }

    emit:
    ch_stlfr2bam
}

process needSplit {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
        
    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path(reads), env(dosplit)

    tag "${id}"

    script:
    def r1 = "${reads[0]}"
    """
    fl=`zcat $r1 | head -1 || true`
    if echo "\$fl" | grep -q '#'; then
        dosplit=0
    else
        dosplit=1
    fi
    """
}
process barcode_split {
    cpus params.CPU1
    memory params.MEM2 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
        
    input:
    tuple val(id), path(reads), val(dosplit)

    when: dosplit

    output:
    tuple val(id), path("${id}_split_1.fq.gz"), path("${id}_split_2.fq.gz"), emit: reads
    tuple val(id), path("split_stat_read1.log"), emit: log

    tag "${id}"
    publishDir "${params.outdir}/$id/fq/", mode: 'link'

    script:
    def bcList = "${params.DB}/barcode/barcode.list"
    def r1 = "${reads[0]}"
    def r2 = "${reads[1]}"
    """
    tmp=`zcat $r1 | head -n 2 | tail -n 1 |wc -c ||true`
    rlen=`expr \$tmp - 1`
    echo \$rlen
    mv $r1 read_1.fq.gz
    mv $r2 read_2.fq.gz

    s1=\$((\$rlen*2+1))
    s2=\$((\$rlen*2+17))
    s3=\$((\$rlen*2+33))

    bcpos="-I \$s1 10 1 false -I \$s2 10 1 false -I \$s3 10 1 false"
    ${params.BIN}MGI.Lite.GenFastQ -F read_1.fq.gz read_2.fq.gz \\
                        -B ${bcList}    \\
                        \$bcpos       \\
                        --stLFR                     \\
                        -O split_out  \\
                        --logPath ./ 

    out1=`find split_out -name "*_1.fq.gz"`
    out2=`find split_out -name "*_2.fq.gz"`
    log=`find split_out -name "split_stat_read.log"`

    mv \$out1 ${id}_split_1.fq.gz
    mv \$out2 ${id}_split_2.fq.gz
    mv \$log split_stat_read1.log
    """
}
process splitRate {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
        
    input:
    tuple val(id), path(reads), val(dosplit)

    when: !dosplit

    output:
    tuple val(id), path("${id}_split_1.fq.gz"), path("${id}_split_2.fq.gz"), emit: reads
    tuple val(id), path("split_stat_read1.log"), emit: log

    tag "${id}"
    publishDir "${params.outdir}/$id/fq/", mode: 'link'

    script:
    def r1 = reads[0]
    def r2 = reads[1]
    """
    mv $r1 ${id}_split_1.fq.gz 2>/dev/null || true
    mv $r2 ${id}_split_2.fq.gz 2>/dev/null || true

    ${params.BIN}python ${params.SCRIPT}/splitRate.py ${id}_split_1.fq.gz > split_stat_read1.log
    """
}
process barcode_split_se600 {
    cpus params.CPU1
    memory params.MEM2 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}_split.fq.gz")

    tag "${id}"
    publishDir "${params.outdir}/$id/fq/", mode: 'link'

    script:
    def bcList = "${params.DB}/barcode/barcode.list"
    def fq = "${reads[0]}"
    """
    ${params.BIN}python3 ${params.SCRIPT}/stlfr2lariatfq_se600.py \\
        $bcList $fq ${id}_split.fq.gz ${params.se600_umi_start}
    """

    stub:
    "touch ${id}_split.fq.gz"
}
process bwa {
    cpus params.cpu3
    memory params.MEM2 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    tag "$id, $lib"

    input:
    val(lib)
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}.${lib}.sort.bam*") 

    // publishDir "${params.outdir}/$id/align/", mode: 'link', enabled: !params.sampleBam
 
    script:
    if (params.ref.startsWith('/')) {
        def ref = params.ref
    } else {
        def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    }
    def r1 = reads[0]
    def r2 = reads[1]
    """
    ${params.BIN}bwa mem -t ${task.cpus} -R '@RG\\tID:${id}\\tSM:sample\\tPL:COMPLETE' $ref $r1 $r2 | \
    ${params.BIN}samtools view -bhS -@ ${task.cpus} -t ${ref}.fai -T $ref - | \
    ${params.BIN}samtools sort -@ ${task.cpus} -T /tmp/sort.${id}.${lib}. -o ${id}.${lib}.sort.bam -
    ${params.BIN}samtools index -@ ${task.cpus} ${id}.${lib}.sort.bam
    """
    stub:
    "touch ${id}.${lib}.sort.bam"
}
process bwaMegabolt {
    cpus params.cpu3
    memory params.MEM2 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.boltq)

    tag "$id, $lib"
    label 'megabolt'

    input:
    val(lib)
    tuple val(id), path(r1), path(r2)

    output:
    tuple val(id), path("${id}.${lib}.*bam*") 

    publishDir "${params.outdir}/$id/align/", mode: 'link', enabled: !params.sampleBam
 
    script:
    if (params.ref.startsWith('/')) {
        def ref = params.ref
        def outbam = "${id}.${lib}.megaboltbwa.bam"
        """
        echo -e \\
        "${id}\\t${r1}\\t${r2}\\t${id}\\t${id}\\t${id}\\tMGISEQ" > ${id}.boltlist

        ${params.MEGABOLT_EXPORT}

        # ${params.MEGABOLT_RUNIT} -l${task.process}.${task.index} ${params.MEGABOLT}                           \\
        ${params.MEGABOLT_RUNIT}  -l\$(basename \$(dirname \$PWD))_\$(basename \$PWD).${task.process}.${task.index} ${params.MEGABOLT} \\
            --type alignmentsortmarkdup --bwa 1  \\
            --ref $ref                              \\
            --list ${id}.boltlist                   \\
            --outputdir .                           

            mv ${id}/${id}.*.bam ${outbam}
            mv ${id}/${id}.*.bam.bai ${outbam}.bai
        """
    } else {
        def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
        def type = "${params.run_bqsr}" ? "alignmentsortmarkdupbqsr" : "alignmentsortmarkdup"
        def gatk = (params.gatk_version == "v4") ? "--gatk4 1" : ""
        def outbam = "${params.run_bqsr}" ? "${id}.${lib}.megaboltbwabqsr.bam" : "${id}.${lib}.megaboltbwa.bam"

        """
        dbsnp=`ls ${params.DB}/${params.ref}/gatk/*dbsnp*.vcf.gz`
        mills=`ls ${params.DB}/${params.ref}/gatk/Mills*.vcf.gz`
        kgsnp=`ls ${params.DB}/${params.ref}/gatk/1000G*snps*.vcf.gz`
        kgindel=`ls ${params.DB}/${params.ref}/gatk/1000G*indels*.vcf.gz`

        echo -e \\
        "${id}\\t${r1}\\t${r2}\\t${id}\\t${id}\\t${id}\\tMGISEQ" > ${id}.boltlist

        ${params.MEGABOLT_EXPORT}

        # ${params.MEGABOLT_RUNIT} -l${task.process}.${task.index} ${params.MEGABOLT}                           \\
        ${params.MEGABOLT_RUNIT}  -l\$(basename \$(dirname \$PWD))_\$(basename \$PWD).${task.process}.${task.index} ${params.MEGABOLT} \\
            --type ${type} --bwa 1  ${gatk}  \\
            --ref $ref                              \\
            --list ${id}.boltlist                   \\
            --vcf \$dbsnp \\
            --knownSites \$dbsnp \\
            --knownSites \$kgindel \\
            --knownSites \$kgsnp \\
            --knownSites \$mills \\
            --outputdir .                           

            mv ${id}/${id}.*.bam ${outbam}
            mv ${id}/${id}.*.bam.bai ${outbam}.bai
        """
    }
    stub:
    "touch ${id}.${lib}.*bam"
}

process kff {    
    cpus params.cpu3
    memory params.MEM2 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    tag "$id"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}.kff")

    // publishDir "${params.outdir}/$id/align/", mode: 'link', enabled: !params.sampleBam

    script:

    """
    printf "%s\\n" ${reads} > file
    kmc -k29 -m${task.memory.giga} -okff -t${task.cpus} @file ${id} .
    """
    stub:
    "touch ${id}.kff"
}
process vg {    
    cpus params.cpu3
    memory params.MEM2 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    tag "$id"

    input:
    tuple val(id), path(kff), path(reads), val(is_pe)

    output:
    tuple val(id), path("${id}.sort.bam*")

    // publishDir "${params.outdir}/$id/align/", mode: 'link', enabled: !params.sampleBam

    script:
    def gbz  = "${params.DB}/hg38/panGenome/hprc-v1.1-mc-grch38.gbz"
    def hapl = "${params.DB}/hg38/panGenome/hprc-v1.1-mc-grch38.hapl"
    def fai  = params.ref.startsWith('/') ? "${params.ref}.fai" : "${params.DB}/hg38/reference/hg38.fa.fai"
    def vg_bin = "/usr/local/app/vg/bin/vg"
    def fq_args = is_pe ? "-f ${reads[0]} -f ${reads[1]}" : "-f ${reads[0]}"
    """
    awk '{print \$1}' $fai | sed 's/^/GRCh38#0#/' > list

    $vg_bin giraffe -Z $gbz --progress --index-basename `pwd`/${id} \\
        --read-group "ID:$id LB:lib1 SM:$id PL:CG PU:unit1" --sample $id -o BAM \\
        --ref-paths list -P -L 3000 $fq_args --kff-name $kff --haplotype-name $hapl \\
        --max-multimaps 3 -t ${task.cpus} | \\
    samtools sort -@ ${task.cpus} -T /tmp/sort.${id}.vg. -o ${id}.sort0.bam -

    samtools view -H ${id}.sort0.bam > header
    sed 's/GRCh38#0#//g' header > new_header.txt
    samtools reheader new_header.txt ${id}.sort0.bam > ${id}.sort.bam

    samtools index -@ ${task.cpus} ${id}.sort.bam
    """
    stub:
    "touch ${id}.sort.bam ${id}.sort.bam.bai"
}

process bqsr {
    cpus params.cpu3
    memory params.MEM2 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    
    input:
    val(lib)
    val(aligner)
    tuple val(id), path(bam) //demo.pf.bwa.merge.bam

    output:
    tuple val(id), path("${id}.${lib}.${aligner}.*.bam")

    tag "$id, $lib, $aligner"
    publishDir "${params.outdir}/$id/align/", mode: 'link'

    script:
    bam = bam.first()
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    cmd = """
    hapmap=${params.DB}/`grep "${params.ref}.hapmap" ${params.DB}/db.list | awk '{print \$2}'`
    dbsnp=${params.DB}/`grep "${params.ref}.dbsnp" ${params.DB}/db.list | awk '{print \$2}'`
    kgindel=${params.DB}/`grep "${params.ref}.1kgindel" ${params.DB}/db.list | awk '{print \$2}'`
    kgsnp=${params.DB}/`grep "${params.ref}.1kgsnp" ${params.DB}/db.list | awk '{print \$2}'`
    mills=${params.DB}/`grep "${params.ref}.1kgmills" ${params.DB}/db.list | awk '{print \$2}'`
    omni=${params.DB}/`grep "${params.ref}.omni" ${params.DB}/db.list | awk '{print \$2}'`
    """
    if (params.gatk_version == "v3") {
        cmd += """
        ${params.BIN}gatk3 -Xmx${task.memory.giga}g -Djava.io.tmpdir=tmpdir \\
        -T BaseRecalibrator \\
        -nct ${task.cpus} \\
        -R $ref \\
        -I $bam \\
        -knownSites \$mills \\
        -knownSites \$kgsnp \\
        -knownSites \$dbsnp \\
        -o sortdup.recal.table

        ${params.BIN}gatk3 -Xmx${task.memory.giga}g -Djava.io.tmpdir=tmpdir \\
        -T PrintReads \\
        -R $ref \\
        -I $bam \\
        -BQSR sortdup.recal.table \\
        -o ${id}.${lib}.${aligner}.bqsr3.bam
        """
    } else if (params.gatk_version == "v4") {
        cmd += """
        ${params.BIN}gatk --java-options "-Xmx${task.memory.giga}g" \\
        BaseRecalibrator \\
        -R $ref \\
        -I $bam \\
        --known-sites \$mills \\
        --known-sites \$kgsnp \\
        --known-sites \$dbsnp \\
        -O sortdup.recal.table

        ${params.BIN}gatk --java-options "-Xmx${task.memory.giga}g" \\
        ApplyBQSR \\
        -R $ref \\
        -I $bam \\
        -bqsr sortdup.recal.table \\
        -O ${id}.${lib}.${aligner}.bqsr4.bam
        """
    }
    return cmd
    stub:
    "touch ${id}.${lib}.${aligner}.bqsr4.bam"
}
process bqsrMegabolt { //stlfr lariat
    label 'megabolt'
    cpus params.cpu3
    memory params.MEM2 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.boltq)

    tag "$id"

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path("${id}.stlfr.lariat.*.bam*") 

    // publishDir "${params.outdir}/$id/align/", mode: 'link'
 
    script:
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    bam = bam.first()
    cmd = """
    kgindel=${params.DB}/`grep "${params.ref}.1kgindel" ${params.DB}/db.list | awk '{print \$2}'`
    kgsnp=${params.DB}/`grep "${params.ref}.1kgsnp" ${params.DB}/db.list | awk '{print \$2}'`
    mills=${params.DB}/`grep "${params.ref}.1kgmills" ${params.DB}/db.list | awk '{print \$2}'`
    dbsnp=${params.DB}/`grep "${params.ref}.dbsnp" ${params.DB}/db.list | awk '{print \$2}'`

    ${params.MEGABOLT_EXPORT}
    """
    if (params.gatk_version == "v3") {
        cmd += """
        ${params.MEGABOLT_RUNIT}  -l\$(basename \$(dirname \$PWD))_\$(basename \$PWD).${task.process}.${task.index} ${params.MEGABOLT} \\
        --type bqsr --bqsr-input $bam  \\
        --ref $ref                              \\
        --vcf \$dbsnp \\
        --knownSites \$dbsnp \\
        --knownSites \$kgindel \\
        --knownSites \$kgsnp \\
        --knownSites \$mills \\
        --outputprefix $id \\
        --outputdir .                           

        mv ${id}/${id}.*.bam ${id}.stlfr.lariat.megaboltbqsr3.bam
        mv ${id}/${id}.*.bai ${id}.stlfr.lariat.megaboltbqsr3.bam.bai                         
        """
    } else if (params.gatk_version == "v4") {
        cmd += """
        ${params.MEGABOLT_RUNIT}  -l\$(basename \$(dirname \$PWD))_\$(basename \$PWD).${task.process}.${task.index} ${params.MEGABOLT} \\
        --type bqsr --bqsr4 1 --bqsr-input $bam  \\
        --ref $ref                              \\
        --vcf \$dbsnp \\
        --knownSites \$dbsnp \\
        --knownSites \$kgindel \\
        --knownSites \$kgsnp \\
        --knownSites \$mills \\
        --outputprefix $id \\
        --outputdir .                           

        mv ${id}/${id}.*.bam ${id}.stlfr.lariat.megaboltbqsr4.bam
        mv ${id}/${id}.*.bai ${id}.stlfr.lariat.megaboltbqsr4.bam.bai                         
        """
    }
    return cmd
    stub:
    "touch ${id}.stlfr.lariat.megaboltbqsr4.bam"
}
process lariatBC {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    
    input:
    tuple val(id), path(r1), path(r2) //hg001_split_1.part_006.fq.gz

    output:
    tuple val(id), path("${id}*.fq.gz")

    tag "$id"
    script:
    cmd = """
    part=`echo $r1 |awk -F '.' '{print \$2}'`
    ${params.BIN}python3 ${params.SCRIPT}/stlfr2lariatfq.py \\
        ${params.DB}/barcode/barcode.list $r1 $r2 ${id}.\$part.fq.gz
    """
    if (!params.keepFiles) {
        cmd += """
        rm `realpath $r1 $r2`
        """
    }
    return cmd
    stub:
    """
    part=`echo $r1 |awk -F '.' '{print \$2}'`
    touch ${id}.\$part.fq.gz
    """
}
process tofake10xHash {
	
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    
	input:
	tuple val(id), val(splitLog)

	output:
	tuple val(id), path("${id}.merge.txt")
	
	tag "$id"
	// publishDir "${params.outdir}/$id/fq/"	//, mode: 'copy'

	script:
	def WL = "${params.DB}/barcode/tenx.whitelist"
	"""
	sed '1,5d' $splitLog | awk '{print \$3,\$2}' > barcode_freq.txt
	perl ${params.SCRIPT}/merge_barcodes.pl barcode_freq.txt  $WL ${id}.merge.txt 1
	"""
    stub:
    "touch ${id}.merge.txt"
}
process tofake10x {
	
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    
	input:
	tuple val(id), path(r1), path(r2), path(hash) //demo2_split_1.9.fq.gz

	output:
	// tuple val(id), path("*map"), emit: stlfr2lariatMap
    tuple val(id), path("*.fake10x_1.fq.gz"), path("*.fake10x_2.fq.gz"), emit: reads
	
	tag "$id"

	script:
	"""
    suffix=`echo $r1 | awk -F "[._]" '{print \$3}'`

	mkdir outdir
	perl ${params.SCRIPT}/fake10x.pl $r1 $r2 $hash outdir \$suffix

    #mv outdir/fake10x_stlfr.map ${r1.getBaseName(2)}.fake10x_stlfr.map
	mv outdir/sample_S1_L001_R1_001.fastq.gz ${r1.getBaseName(2)}.fake10x_1.fq.gz
	mv outdir/sample_S1_L001_R2_001.fastq.gz ${r2.getBaseName(2)}.fake10x_2.fq.gz
	"""
    stub:
    "touch ${r1.getBaseName(2)}.fake10x_1.fq.gz ${r1.getBaseName(2)}.fake10x_2.fq.gz"
}
process mergeMaps {
	
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    
	input:
	tuple val(id), path(maps)

	output:
	tuple val(id), path("${id}.fake10xmap"), emit: stlfr2lariatMap
	
	tag "$id"
	// publishDir "${params.outdir}/$id/fq/"	//, mode: 'copy'

	script:
	"""
	cat $maps > ${id}.fake10xmap
	"""
}

process fake10x2lariat {
    
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(fq1), path(fq2) //demo2.stlfr.005_1.fake10x_1.fq.gz

    output:
    tuple val(id), path("*.lariat.fq.gz") //demo2.stlfr.005_1.lariat.fq.gz

    tag "$id"
 
    script:
    """
    python3 ${params.SCRIPT}/stlfr2lariat_v3.py $fq1 $fq2 ${fq1.getBaseName(3)}.lariat.fq.gz 
    """
    stub:
    "touch ${fq1.getBaseName(3)}.lariat.fq.gz "
}
process mergeFq {
    
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(fqs) //demo2.stlfr.005_1.lariat.fq.gz

    output:
    tuple val(id), path("${id}.lariat.fq.gz")

    tag "$id"
    // publishDir "${params.outdir}/$id/fq/"
 
    script:
    cmd = """
	echo `ls *.fq.gz| sort -t '.' -k 2n` > tmp
	cat `ls *.fq.gz| sort -t '.' -k 2n`  > ${id}.lariat.fq.gz
    # ${params.BIN}seqkit scat -j ${task.cpus} . -g 
    """
    if (!params.keepFiles) {
        cmd += """
        find . -type l -name "*.fq.gz" -exec sh -c 'rm "\$(readlink -f "{}")"' \\;
        """
    }
    return cmd
    stub:
    "touch ${id}.lariat.fq.gz"
}

process lariat {  
    cpus params.cpu2
    memory params.MEM3 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(fq)

    output:
    tuple val(id), path("bc_sorted_bam.bam")
    tag "$id"
    // publishDir "${params.outdir}/$id/align/"
 
    script:
    def ref = params.ref.startsWith('/') ? params.ref : "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    def t = "${params.lariatStLFRBC}" ? 0 : 7
    cmd = """
    ${params.BIN}lariat -genome $ref -output . -reads $fq -threads=${task.cpus} -trim_length $t -read_groups "1:N:0:NAAGTGCT:0"
    # ${params.BIN}samtools sort -@ ${task.cpus} -m 1G -o bc_re-sorted_bam.bam bc_sorted_bam.bam

    rm -f 000* ZZZ*

    # mv bc_re-sorted_bam.bam ${id}.lariat.sort.bam 
    # ${params.BIN}samtools index ${id}.lariat.sort.bam
    """
    if (!params.keepFiles) {
        cmd += """
        rm -f `realpath $fq`
        """
    }
    return cmd
    stub:
    "touch bc_sorted_bam.bam"
}
process sortbam {  
    cpus params.cpu2
    memory params.MEM3 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path("${id}.lariat.sort.bam*")
    tag "$id"
    // publishDir "${params.outdir}/$id/align/"
 
    script:
    def ref = params.ref.startsWith('/') ? params.ref : "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    def t = "${params.lariatStLFRBC}" ? 0 : 7
    cmd = """
    ${params.BIN}samtools sort -@ ${task.cpus} -m 1G -o bc_re-sorted_bam.bam $bam

    mv bc_re-sorted_bam.bam ${id}.lariat.sort.bam 
    ${params.BIN}samtools index ${id}.lariat.sort.bam
    """
    if (!params.keepFiles) {
        cmd += """
        rm -f `realpath $bam`
        """
    }
    return cmd
    stub:
    "touch ${id}.lariat.sort.bam"
}
process markdup {
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    val(lib)
    val(aligner)
    tuple val(id), path(idxbam)

    output:
    tuple val(id), path("${id}.${lib}.${aligner}.*.bam*")

    tag "$id"

    publishDir "${params.outdir}/$id/align/", mode: 'link', enabled: !params.sampleBam

    script:
    bam = idxbam.first()
    cmd = ""
    if (params.markdup == "sambamba") {
        cmd += """
        ${params.BIN}sambamba markdup \\
            -t ${task.cpus} --tmpdir /tmp $bam ${id}.${lib}.${aligner}.sambamba.bam
        ${params.BIN}samtools index -@ ${task.cpus} ${id}.${lib}.${aligner}.sambamba.bam
        """
    } else if (params.markdup == "picard") {
        cmd += """
        java -Xms${task.memory.giga}g -Xmx${task.memory.giga}g -jar ${params.SCRIPT}/picard/picard.jar MarkDuplicates I=$bam O=${id}.${lib}.${aligner}.picard.bam M=${id}.${lib}.${aligner}.picardMarkdup.log TMP_DIR=/tmp
        ${params.BIN}samtools index -@ ${task.cpus} ${id}.${lib}.${aligner}.picard.bam
        """
    } else if (params.markdup == "biobambam2") {
        cmd += """
        ${params.BIN}bammarkduplicates2 I=$bam O=${id}.${lib}.${aligner}.biobambam2.bam M=${id}.${lib}.${aligner}.biobambam2.log markthreads=${task.cpus}
        ${params.BIN}samtools index -@ ${task.cpus} ${id}.${lib}.${aligner}.biobambam2.bam
        """
    } else if (params.markdup == "gatk4") {
        cmd += """
        ${params.BIN}gatk MarkDuplicatesSpark \\
            --spark-master local[${task.cpus}] \\
            -I $bam -O ${id}.${lib}.${aligner}.MarkDuplicatesSpark.bam -M ${id}.${lib}.${aligner}.MarkDuplicatesSpark.log

        ${params.BIN}samtools index -@ ${task.cpus} ${id}.${lib}.${aligner}.MarkDuplicatesSpark.bam
        """
    }
    if (!params.keepFiles) {
        cmd += """
        rm `realpath $idxbam`
        """
    }
    return cmd
    stub:
    "touch ${id}.${lib}.${aligner}.${params.markdup}.bam ${id}.${lib}.${aligner}.${params.markdup}.bam.bai"
}
process sampleBam_samtools {
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    val(lib)
    val(aligner)
    tuple val(id), path(indexedbam)

    output:
    tuple val(id), path("${id}.${lib}.${aligner}.sampled.bam*"), emit: bam
    //path("${id}.${lib}.${aligner}.sampledbamdepth")

    tag "$id, $lib, $aligner"
    publishDir "${params.outdir}/$id/align/", mode: 'link'
 
    script:
    def bam = indexedbam.first()
    def cov = lib == "stlfr" ? "${params.stLFR_sampling_cov}" : "${params.PF_sampling_cov}"
    if (params.debug) {
        cmd = """
        cp $bam ${id}.${lib}.${aligner}.sampled.bam
        cp ${bam}.bai ${id}.${lib}.${aligner}.sampled.bam.bai
        """
    } else {
        cmd = """
        bamcov=`${params.BIN}samtools depth -@ ${task.cpus} $bam | awk '{sum += \$3}END{print sum/${params.ref_len}}'`
        echo \$bamcov > tmp
        ratio=`echo "scale=5; $cov/\$bamcov" | bc`
        if [[ \$ratio > 1 ]];then
            ratio=1.0
        fi

        ${params.BIN}samtools view -@ ${task.cpus} -s \$ratio $bam -o ${id}.${lib}.${aligner}.sampled.bam
        ${params.BIN}samtools index -@ ${task.cpus} ${id}.${lib}.${aligner}.sampled.bam
        """
    }

    if (!params.keepFiles) {
        cmd += """
        rm `realpath $indexedbam`
        """
    }
    return cmd
}
process sampleBam {
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    val(lib)
    val(aligner)
    tuple val(id), path(indexedbam)

    output:
    tuple val(id), path("${id}.${lib}.${aligner}.sampled.bam*"), emit: bam
    //path("${id}.${lib}.${aligner}.sampledbamdepth")

    tag "$id, $lib, $aligner"
    publishDir "${params.outdir}/$id/align/", mode: 'link'
 
    script:
    def bam = indexedbam.first()
    def cov = lib == "stlfr" ? "${params.stLFR_sampling_cov}" : "${params.PF_sampling_cov}"

    cmd = """
    bamcov=`${params.BIN}samtools depth -@ ${task.cpus} $bam | awk '{sum += \$3}END{print sum/${params.ref_len}}'`
    echo \$bamcov > tmp
    ratio=`echo "scale=5; $cov/\$bamcov" | bc`
    if [[ \$ratio > 1 ]];then
        cp $bam ${id}.${lib}.${aligner}.sampled.bam
        cp ${bam}.bai ${id}.${lib}.${aligner}.sampled.bam.bai
    else 
        bamcov=`samtools stats $bam | awk -v ref=${params.ref_len} '\$2=="bases" && \$3=="mapped" && \$4=="(cigar):" {print \$5/ref}'`
        echo \$bamcov > tmp
        ratio=`echo "scale=5; $cov/\$bamcov" | bc`
        if [[ \$ratio > 1 ]];then
            cp $bam ${id}.${lib}.${aligner}.sampled.bam
            cp ${bam}.bai ${id}.${lib}.${aligner}.sampled.bam.bai
        else
            seed=42
            ${params.BIN}samtools view -@ ${task.cpus} -s \$seed\$ratio -b $bam > ${id}.${lib}.${aligner}.sampled.bam
            ${params.BIN}samtools index -@ ${task.cpus} ${id}.${lib}.${aligner}.sampled.bam
        fi
    fi
    """
    if (!params.keepFiles) {
        cmd += """
        rm `realpath $indexedbam`
        """
    }
    return cmd
    stub:
    "touch ${id}.${lib}.${aligner}.sampled.bam"
}

process depth {
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    val(lib)
    val(aligner)
    tuple val(id), path(bam)

    output:
    tuple val(id), path("*.bed") //demo.stlfr.bwa.cov10.bed

    tag "$id, $lib, $aligner"
    // publishDir "${params.outdir}/$id/align/"
    
    script:
    def bam = bam.first() //demo.stlfr.bwa.bam
    s = lib == "stlfr" ? ">" : "<"
    """
    ${params.BIN}samtools depth -@ ${task.cpus} $bam | awk -v cov="${params.bamcov}" '\$3 $s cov {print \$1"\\t"\$2"\\t"\$2}' | \\
    ${params.BIN}bedtools merge -i > ${id}.${lib}.${aligner}.cov${params.bamcov}.bed

    """
    stub:
    "touch ${id}.${lib}.${aligner}.cov${params.bamcov}.bed"
}

process bed {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    val(aligner)
    tuple val(id), path(bed1), path(bed2) //demo.stlfr.bwa.cov10.bed

    output:
    tuple val(id), path("*intersect.bed") 

    tag "$id, $aligner"
    // publishDir "${params.outdir}/$id/align/", mode: 'link'
 
    script:
    """
    ${params.BIN}bedtools intersect -a $bed1 -b $bed2 > ${id}.${aligner}.cov${params.bamcov}.intersect.bed
    """
    stub:
    "touch ${id}.${aligner}.cov${params.bamcov}.intersect.bed"
}



process addBxSe600 {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(bam_files)

    output:
    tuple val(id), path("${id}.stlfr2.vg.bx.bam*")

    tag "$id"

    script:
    def bam = bam_files instanceof List ? bam_files.first() : bam_files
    """
    python3 ${params.SCRIPT}/add_bx_se600.py $bam ${id}.stlfr2.vg.bx.bam
    ${params.BIN}samtools index -@ ${task.cpus} ${id}.stlfr2.vg.bx.bam
    """
    stub:
    "touch ${id}.stlfr2.vg.bx.bam ${id}.stlfr2.vg.bx.bam.bai"
}

process stLFRQC {
    cpus params.CPU1
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path("${id}.lfr.report"), emit: report
    path("06*.txt")

    tag "$id"
    publishDir "${params.outdir}/report/$id/", mode:'link'
 
    script:
    bam = bam.first()
    def ref = params.ref.startsWith('/') ? params.ref : "${params.DB}/${params.ref}/reference/${params.ref}.fa"    
    """
    ${params.BIN}stLFRQC --samtools /usr/bin/samtools --python3 /usr/bin/python3 \
        --ref $ref \
        --bam $bam \
        --thread ${task.cpus} \
        --script ${params.SCRIPT}/stLFRQC || echo "stLFRQC failed, using NA placeholders" >&2

    if [ -f 06.lfr_highquality.txt ]; then
        sed -n '11p' 06.lfr_highquality.txt > tmp
        sed -n '4p' 06.lfr_length.txt >> tmp
        sed -n '4p' 06.lfr_readpair.txt >> tmp
        sed -n '77p' 06.lfr_per_barcode.txt >> tmp
    else
        printf 'NA\nNA\nNA\nNA\n' > tmp
    fi
    mv tmp ${id}.lfr.report
    """
    stub:
    "touch 06*.txt ${id}.lfr.report "
}
