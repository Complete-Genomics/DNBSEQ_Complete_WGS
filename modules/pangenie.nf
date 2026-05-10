process pangenie {
    cpus params.cpu3
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    when: params.ref == 'hg38' || params.ref.contains('GRCh38')

    input:
    tuple val(id), path(r)

    output:
    tuple val(id), path("pangenie_genotyping_biallelic.vcf.gz*")

    tag "$id"
    publishDir "${params.outdir}/$id/", mode: 'link'

    script:
    def python = "/usr/local/app/miniconda3/bin/python"
    def py = "/usr/local/app/pangenie/pipelines/run-from-callset/scripts/convert-to-biallelic.py"
    """
    cat ${r[1]} ${r[2]} | gunzip > merge.fq
    PanGenie -f ${params.DB}/pangenie/HPRC_index -i merge.fq -o pangenie -j ${task.cpus} -t ${task.cpus}

    cat pangenie_genotyping.vcf | $python $py ${params.DB}/pangenie/cactus_filtered_ids_biallelic.vcf.gz |bgzip > pangenie_genotyping_biallelic.vcf.gz
    tabix pangenie_genotyping_biallelic.vcf.gz

    rm merge.fq pangenie_genotyping.vcf
    """
}
process pangenie_plot {
    cpus params.cpu3
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(vcf), path(hapblock)

    output:
    tuple val(id), path("chromosome_sv.png")

    // cache false
    tag "$id"
    publishDir "${params.outdir}/report/$id/"

    script:
    vcf = vcf.first()
    """
    bcftools view -H \
    -e 'GT="0/0" || GT="./." || GT="./0" || GT="0/." || GT="." || GT="0"' \
    $vcf |
    awk -F'\\t' '
    BEGIN{OFS="\\t"; print "Type","Shape","Chr","Start","End","color"}
    {
        split(\$3,a,/>/); L=a[4]; R=a[5]; len=R-L
        if(len<=10000) next

        ref=length(\$4); alt=length(\$5)
        if(\$5=="<DEL>" || (ref>1 && alt==1))        t="Deletion"
        else if(\$5=="<INS>" || (ref==1 && alt>1))  t="Insertion"
        else                                       t="Complex"

        shape = (t=="Complex"?"box":(t=="Deletion"?"triangle":"circle"))
        color = (t=="Complex"?"6a3d9a":(t=="Deletion"?"ff7f01":"33a02c"))
        chr=\$1; sub(/^chr/,"", chr)
        print t, shape, chr, \$2, \$2+1, color
    }' > sv_10k.txt

    python ${params.SCRIPT}/band.py $hapblock
    Rscript ${params.SCRIPT}/pangenie_plot.R sv_10k.txt
    convert -crop 100x66%+0+0 chromosome.png chromosome_sv.png
    """
}
process pangenie_var_plot {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(vcf)

    output:
    path "pangenie_var_plot.png"

    tag "$id"
    publishDir "${params.outdir}/report/$id/", mode: 'link'

    script:
    """
    python ${params.SCRIPT}/pangenie_var_plot.py $vcf
    """
}
process pangenie_frombam {
    cpus params.cpu3
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    when: params.ref == 'hg38' || params.ref.contains('GRCh38')

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path("pangenie_genotyping_biallelic.vcf.gz*")

    tag "$id"
    publishDir "${params.outdir}/$id/", mode: 'link'

    script:
    def bam_file = bam.first()
    def python = "/usr/local/app/miniconda3/bin/python"
    def py = "/usr/local/app/pangenie/pipelines/run-from-callset/scripts/convert-to-biallelic.py"
    """
    ${params.BIN}samtools fastq -@ ${task.cpus} -0 /dev/null -s /dev/null \
        -1 reads_r1.fq.gz -2 reads_r2.fq.gz $bam_file

    cat reads_r1.fq.gz reads_r2.fq.gz | gunzip > merge.fq

    PanGenie -f ${params.DB}/pangenie/HPRC_index -i merge.fq -o pangenie -j ${task.cpus} -t ${task.cpus}

    cat pangenie_genotyping.vcf | $python $py ${params.DB}/pangenie/cactus_filtered_ids_biallelic.vcf.gz | bgzip > pangenie_genotyping_biallelic.vcf.gz
    tabix pangenie_genotyping_biallelic.vcf.gz

    rm merge.fq reads_r1.fq.gz reads_r2.fq.gz pangenie_genotyping.vcf
    """
}
