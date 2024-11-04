process bwa {    
    cpus params.cpu3
    memory params.MEM2 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

    tag "$id, $lib"

    input:
    val(lib)
    tuple val(id), path(r1), path(r2)

    output:
    tuple val(id), path("${id}.${lib}.sort.bam*") 

    // publishDir "${params.outdir}/$id/align/", mode: 'link', enabled: !params.sampleBam
 
    script:
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    """
    ${params.BIN}bwa mem -t ${task.cpus} -R '@RG\\tID:${id}\\tSM:sample\\tPL:COMPLETE' $ref $r1 $r2 | \
    ${params.BIN}samtools view -bhS -@ ${task.cpus} -t ${ref}.fai -T $ref - | \
    ${params.BIN}samtools sort -@ ${task.cpus} -T `pwd`/sort.tmp. -o ${id}.${lib}.sort.bam - 
    ${params.BIN}samtools index -@ ${task.cpus} ${id}.${lib}.sort.bam
    """
    stub:
    "touch ${id}.${lib}.sort.bam"
}
process bwaMegabolt {
    cpus params.cpu3
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.boltq} ${params.extraCluOpt}"

    tag "$id, $lib"
    label 'megabolt'

    input:
    val(lib)
    tuple val(id), path(r1), path(r2)

    output:
    tuple val(id), path("${id}.${lib}.*bam*") 

    publishDir "${params.outdir}/$id/align/", mode: 'link', enabled: !params.sampleBam
 
    script:
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    def refdir = "${params.DB}/${params.ref}/reference/"
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

    ${params.MEGABOLT_RUNIT} -l${task.process}.${task.index} ${params.MEGABOLT}                           \\
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
    stub:
    "touch ${id}.${lib}.*bam"
}

process bqsr {
    cpus params.cpu3
    memory params.MEM2 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
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
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.boltq} ${params.extraCluOpt}"

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
        ${params.MEGABOLT_RUNIT} -l${task.process}.${task.index} ${params.MEGABOLT}                           \\
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
        ${params.MEGABOLT_RUNIT} -l${task.process}.${task.index} ${params.MEGABOLT}                           \\
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
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
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
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
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
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
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
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
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
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

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
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

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
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

    input:
    tuple val(id), path(fq)

    output:
    tuple val(id), path("bc_sorted_bam.bam")
    tag "$id"
    // publishDir "${params.outdir}/$id/align/"
 
    script:
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
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
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path("${id}.lariat.sort.bam*")
    tag "$id"
    // publishDir "${params.outdir}/$id/align/"
 
    script:
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
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
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

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
            -t ${task.cpus} --tmpdir . $bam ${id}.${lib}.${aligner}.sambamba.bam
        ${params.BIN}samtools index -@ ${task.cpus} ${id}.${lib}.${aligner}.sambamba.bam
        """
    } else if (params.markdup == "picard") {
        cmd += """
        java -Xms${task.memory.giga}g -Xmx${task.memory.giga}g -jar ${params.SCRIPT}/picard/picard.jar MarkDuplicates I=$bam O=${id}.${lib}.${aligner}.picard.bam M=${id}.${lib}.${aligner}.picardMarkdup.log TMP_DIR=. 
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
    "touch ${id}.${lib}.${aligner}.*.bam"
}
process sampleBam_samtools {
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

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
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

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
        ${params.BIN}gatk DownsampleSam -I $bam -O ${id}.${lib}.${aligner}.sampled.bam -P \$ratio
        ${params.BIN}samtools index -@ ${task.cpus} ${id}.${lib}.${aligner}.sampled.bam
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
process intersect {
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

    input:
    val(aligner)
    tuple val(id), path(stlfrbam), path(pfbam)

    output:
    tuple val(id), path("*.intersect.bed") //demo.stlfr.bwa.cov10.bed

    tag "$id, $aligner"
    publishDir "${params.outdir}/$id/align/", mode: 'link'
    
    script:
    stlfrbam = stlfrbam.first() //demo.stlfr.bwa.bam
    pfbam = pfbam.first()
    """
    ${params.BIN}samtools depth -@ ${task.cpus} $stlfrbam $pfbam | awk -v cov="${params.PF_lt_stLFR_depth}" '\$3 >= cov && \$4 < cov {print \$1"\\t"\$2"\\t"\$2}' | \\
    ${params.BIN}bedtools merge > ${id}.${aligner}.cov${params.PF_lt_stLFR_depth}.intersect.bed

    """
    stub:
    "touch ${id}.${aligner}.cov${params.PF_lt_stLFR_depth}.intersect.bed"
}
process depth {
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

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
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

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

process mergeBam {
    cpus params.cpu2
    memory params.MEM3 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

    input:
    val(aligner)
    tuple val(id), path(pfbam), path(stlfrbam), path(bed) 

    output:
    tuple val(id), path("${id}.${aligner}.merge.bam*") 
    tag "$id, $aligner"

    publishDir "${params.outdir}/$id/align/", mode: 'link'

    script:
    def pfbam = pfbam.first()
    def stlfrbam = stlfrbam.first()
    """
    ${params.BIN}samtools view -hb -L $bed $stlfrbam > lfr_lfr10_pf10.bam && \\
    ${params.BIN}samtools index lfr_lfr10_pf10.bam && \\
    ${params.BIN}samtools addreplacerg -w -O BAM -r '@RG\\tID:${id}\\tSM:sample' -o new.bam lfr_lfr10_pf10.bam && \\
    ${params.BIN}samtools index new.bam && \\
    ${params.BIN}samtools reheader $pfbam new.bam > new2.bam && \\
    ${params.BIN}samtools index new2.bam && \\
    # ${params.BIN}samtools addreplacerg -O BAM -r '@RG\\tID:sample\\tSM:sample' -o new3.bam new2.bam && \\
    # ${params.BIN}samtools index new3.bam && \\
    ${params.BIN}samtools merge -@ ${task.cpus} -f ${id}.${aligner}.merge.bam $pfbam new2.bam && \\
    ${params.BIN}samtools index -@ ${task.cpus} ${id}.${aligner}.merge.bam

    """
    stub:
    "touch ${id}.${aligner}.merge.bam"
}

process stLFRQC {
    cpus params.CPU1
    memory params.MEM1 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path("${id}.lfr.report"), emit: report
    path("06*.txt")

    tag "$id"
    publishDir "${params.outdir}/report/$id/", mode:'link'
 
    script:
    bam = bam.first()
    def ref = "${params.DB}/${params.ref}/reference/${params.ref}.fa"    
    """
    ${params.BIN}stLFRQC --samtools /usr/bin/samtools --python3 /usr/bin/python3 \
        --ref $ref \
        --bam $bam \
        --thread ${task.cpus}

    sed -n '11p' 06.lfr_highquality.txt > tmp
    sed -n '4p' 06.lfr_length.txt >> tmp
    sed -n '4p' 06.lfr_readpair.txt >> tmp
    sed -n '77p' 06.lfr_per_barcode.txt >> tmp
    mv tmp  ${id}.lfr.report
    """
    stub:
    "touch 06*.txt ${id}.lfr.report "
}