// per-haplotype overlapping-window de novo assembly for stLFR
//
// input:  [id, merged.bam, phased.vcf.gz]
// output: per sample, hp1/hp2 contigs.fa.gz
//
// requires (must be on PATH inside container or BIN):
//   whatshap, samtools, bedtools, megahit, bgzip


process haplotag {
    cpus params.cpu3
    memory params.MEM2 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(bam), path(vcf)

    output:
    tuple val(id), path("${id}.haplotag.bam"), path("${id}.haplotag.bam.bai"), emit: bam

    tag "$id"
    publishDir "${params.outdir}/$id/haplodenovo/", mode: 'link'

    script:
    def ref      = "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    def vcf_file = vcf instanceof List ? vcf.first() : vcf
    def bam_file = bam instanceof List ? bam.first() : bam
    """
    ${params.BIN}whatshap haplotag \\
        --reference $ref \\
        --output ${id}.haplotag.bam \\
        --ignore-read-groups \\
        --skip-missing-contigs \\
        --output-threads ${task.cpus} \\
        $vcf_file $bam_file

    ${params.BIN}samtools index -@ ${task.cpus} ${id}.haplotag.bam
    """
    stub:
    "touch ${id}.haplotag.bam ${id}.haplotag.bam.bai"
}


process propagateBxHp {
    cpus params.cpu3
    memory params.MEM2 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(bam), path(bai)

    output:
    tuple val(id), path("${id}.bcprop.bam"), path("${id}.bcprop.bam.bai"), emit: bam
    path("${id}.bx_hp.tsv")

    tag "$id"
    publishDir "${params.outdir}/$id/haplodenovo/", mode: 'link'

    script:
    // ratio: BX must have majority HP at >=5x over the other HP to be assigned
    """
    # pass1: BX majority vote -> bx_hp.tsv
    ${params.BIN}samtools view -@ ${task.cpus} ${bam} \\
    | awk 'BEGIN{OFS="\\t"}
        {
            bx=""; hp=""
            for(i=12;i<=NF;i++){
                if(substr(\$i,1,5)=="BX:Z:") bx=substr(\$i,6)
                else if(substr(\$i,1,5)=="HP:i:") hp=substr(\$i,6)
            }
            if(bx!="" && (hp=="1" || hp=="2")) cnt[bx,hp]++
        }
        END{
            for(k in cnt){
                split(k,a,SUBSEP); b=a[1]; h=a[2]
                if(seen[b]) continue
                seen[b]=1
                c1=cnt[b,"1"]+0; c2=cnt[b,"2"]+0
                if(c1 >= 5*c2 && c1 > 0) print b, "1"
                else if(c2 >= 5*c1 && c2 > 0) print b, "2"
            }
        }' | sort -u > ${id}.bx_hp.tsv

    # pass2: re-tag all reads using BX -> HP map
    ${params.BIN}samtools view -h -@ ${task.cpus} ${bam} \\
    | python3 ${params.SCRIPT}/propagate_bx_hp.py --map ${id}.bx_hp.tsv \\
    | ${params.BIN}samtools view -bS -@ ${task.cpus} - > ${id}.bcprop.bam

    ${params.BIN}samtools index -@ ${task.cpus} ${id}.bcprop.bam
    """
    stub:
    "touch ${id}.bcprop.bam ${id}.bcprop.bam.bai ${id}.bx_hp.tsv"
}


process makeWindows {
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    val(id)

    output:
    path("windows.bed")

    tag "$id"

    script:
    def fai  = "${params.DB}/${params.ref}/reference/${params.ref}.fa.fai"
    def win  = params.denovo_window_size ?: 60000
    def step = params.denovo_window_step ?: 30000
    """
    awk 'BEGIN{OFS="\\t"} \$1 ~ /^chr([0-9]+|X|Y)\$/ {print \$1, 0, \$2}' $fai \\
    | ${params.BIN}bedtools makewindows -b - -w ${win} -s ${step} > windows.bed
    """
    stub:
    "touch windows.bed"
}


process denovoBatch {
    cpus params.cpu2
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    errorStrategy 'ignore'

    input:
    tuple val(id), path(bam), path(bai), val(hp), path(bed_chunk)

    output:
    tuple val(id), val(hp), path("contigs.${bed_chunk.name}.hp${hp}.fa.gz"), optional: true

    tag "$id hp${hp} ${bed_chunk.name}"

    script:
    def min_reads   = params.denovo_min_reads   ?: 50
    def min_contig  = params.denovo_min_contig  ?: 500
    def max_reads   = params.denovo_max_reads   ?: 20000  // subsample cap per window
    """
    > out.fa
    while IFS=\$'\\t' read -r chr s e; do
        s1=\$((s+1))
        region="\${chr}:\${s1}-\${e}"
        safe="\${chr}_\${s}_\${e}"

        ${params.BIN}samtools view -b -d HP:${hp} -@ ${task.cpus} ${bam} \$region > w.bam 2>/dev/null
        n=\$(${params.BIN}samtools view -c w.bam)
        if [ "\$n" -lt ${min_reads} ]; then rm -f w.bam; continue; fi

        # subsample if over cap
        if [ "\$n" -gt ${max_reads} ]; then
            frac=\$(awk -v n=\$n -v c=${max_reads} 'BEGIN{printf "%.4f", c/n}')
            ${params.BIN}samtools view -b -s 42\${frac#0} -@ ${task.cpus} w.bam > w_s.bam
            mv w_s.bam w.bam
        fi

        ${params.BIN}samtools fastq -@ ${task.cpus} \\
            -1 r1.fq -2 r2.fq -0 /dev/null -s /dev/null w.bam 2>/dev/null

        rm -rf asm
        megahit -1 r1.fq -2 r2.fq -o asm -t ${task.cpus} \\
            --min-contig-len ${min_contig} --memory 0.2 -f >/dev/null 2>&1 || true

        if [ -s asm/final.contigs.fa ]; then
            awk -v r=\$safe -v h=${hp} \\
                '/^>/{print ">"r".hp"h".c"++c; next}{print}' \\
                asm/final.contigs.fa >> out.fa
        fi
        rm -f r1.fq r2.fq w.bam
        rm -rf asm
    done < ${bed_chunk}

    if [ -s out.fa ]; then
        ${params.BIN}bgzip -@ ${task.cpus} -c out.fa > contigs.${bed_chunk.name}.hp${hp}.fa.gz
    fi
    rm -f out.fa
    """
    stub:
    "touch contigs.${bed_chunk.name}.hp${hp}.fa.gz"
}


process catContigs {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), val(hp), path(fa_gzs)

    output:
    tuple val(id), val(hp), path("${id}.hp${hp}.contigs.fa.gz"), emit: fa
    path("${id}.hp${hp}.contigs.stats.txt")

    tag "$id hp${hp}"
    publishDir "${params.outdir}/$id/haplodenovo/", mode: 'link'

    script:
    """
    cat ${fa_gzs} > ${id}.hp${hp}.contigs.fa.gz

    # quick N50 / count
    ${params.BIN}samtools faidx ${id}.hp${hp}.contigs.fa.gz 2>/dev/null || \\
        zcat ${id}.hp${hp}.contigs.fa.gz | awk '/^>/{n++; next}{l+=length(\$0)} END{print "contigs="n"\\ttotal_bp="l}' \\
            > ${id}.hp${hp}.contigs.stats.txt

    if [ ! -s ${id}.hp${hp}.contigs.stats.txt ]; then
        zcat ${id}.hp${hp}.contigs.fa.gz \\
        | awk '/^>/{if(s)print length(s); s=""; n++; next}{s=s\$0} END{if(s)print length(s); print "n_contigs="n > "/dev/stderr"}' \\
        | sort -nr > lens.txt
        total=\$(awk '{s+=\$1}END{print s}' lens.txt)
        awk -v T=\$total 'BEGIN{half=T/2}{c+=\$1; if(c>=half){print "N50="\$1"\\ttotal_bp="T; exit}}' lens.txt > ${id}.hp${hp}.contigs.stats.txt
        rm -f lens.txt
    fi
    """
    stub:
    "touch ${id}.hp${hp}.contigs.fa.gz ${id}.hp${hp}.contigs.stats.txt"
}


workflow WF_haplodenovo {
    take:
    ch_in   // [id, bam, phased_vcf]

    main:
    haplotag(ch_in)
    propagateBxHp(haplotag.out.bam).set { ch_prop }

    // single windows.bed for the run (ref-based, sample-independent)
    makeWindows(ch_prop.bam.map{ it[0] }.first()).set { ch_winbed }

    // split BED into chunks
    def chunk = params.denovo_chunk_lines ?: 200
    ch_winbed
        .splitText(by: chunk, file: true)
        .set { ch_chunks }

    // [id, bam, bai, hp, chunk]
    ch_prop.bam
        .combine( Channel.fromList(['1','2']) )
        .combine( ch_chunks )
        .set { ch_jobs }

    denovoBatch(ch_jobs).set { ch_batches }

    ch_batches.groupTuple(by: [0,1]).set { ch_grouped }
    catContigs(ch_grouped).set { ch_final }

    emit:
    contigs = ch_final.fa
}