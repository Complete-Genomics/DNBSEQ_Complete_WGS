process pangenie {
    cpus params.cpu3
    memory params.pangenie_memory
    maxForks 2
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    when: params.ref == 'hg38' || params.ref.contains('GRCh38')

    input:
    tuple val(id), path(r)

    output:
    tuple val(id), path("{pangenie_genotyping_biallelic.vcf.gz*,pangenie.FAILED}")

    tag "$id"
    publishDir "${params.outdir}/$id/", mode: 'link'

    script:
    def python = "/usr/local/app/miniconda3/bin/python"
    def py = "/usr/local/app/pangenie/pipelines/run-from-callset/scripts/convert-to-biallelic.py"
    """
    # Never fail the task: the final report must still be generated.
    # On failure write pangenie.FAILED + empty placeholder VCF; downstream plots then emit *.FAILED.
    set +e
    (
    set -e
        cat ${r.join(' ')} | gunzip > merge.fq
        PanGenie -f ${params.DB}/pangenie/HPRC_index -i merge.fq -o pangenie -j ${task.cpus} -t ${task.cpus}

        cat pangenie_genotyping.vcf | $python $py ${params.DB}/pangenie/cactus_filtered_ids_biallelic.vcf.gz |bgzip > pangenie_genotyping_biallelic.vcf.gz
        tabix pangenie_genotyping_biallelic.vcf.gz
    )
    rc=\$?
    set -e
    rm -f merge.fq pangenie_genotyping.vcf reads_r1.fq.gz reads_r2.fq.gz
    if [ \$rc -ne 0 ]; then
        echo "PanGenie exited with status \$rc" | tee pangenie.FAILED >&2
        : > pangenie_genotyping_biallelic.vcf.gz
        : > pangenie_genotyping_biallelic.vcf.gz.tbi
    fi
    exit 0
    """
    stub:
    "touch pangenie_genotyping_biallelic.vcf.gz pangenie_genotyping_biallelic.vcf.gz.tbi"
}
process pangenie_plot {
    cpus params.cpu3
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(vcf), path(hapblock)

    output:
    tuple val(id), path("chromosome_sv.*")   // chromosome_sv.png, or chromosome_sv.FAILED

    // cache false
    tag "$id"
    publishDir "${params.outdir}/report/$id/", mode: 'copy'

    stub:
    "touch chromosome_sv.png"

    script:
    vcf = vcf.first()
    """
    # Soft-fail: write chromosome_sv.FAILED (read by scripts/my_html.py) instead of failing the task.
    if [ -e pangenie.FAILED ] || [ ! -s $vcf ]; then
        { cat pangenie.FAILED 2>/dev/null || echo "upstream PanGenie output missing or empty"; } > chromosome_sv.FAILED
        exit 0
    fi
    set +e
    (
    set -e
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
    )
    rc=\$?
    set -e
    if [ \$rc -ne 0 ] || [ ! -s chromosome_sv.png ]; then
        echo "pangenie_plot exited with status \$rc" > chromosome_sv.FAILED
        rm -f chromosome_sv.png
    fi
    exit 0
    """
}
process pangenie_var_plot {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(vcf)

    output:
    path "pangenie_var_plot.*"   // pangenie_var_plot.png, or pangenie_var_plot.FAILED

    tag "$id"
    publishDir "${params.outdir}/report/$id/", mode: 'copy'

    stub:
    "touch pangenie_var_plot.png"

    script:
    def vcf0 = vcf instanceof List ? vcf.first() : vcf
    """
    # Soft-fail: write pangenie_var_plot.FAILED (read by scripts/my_html.py) instead of failing the task.
    if [ -e pangenie.FAILED ] || [ ! -s $vcf0 ]; then
        { cat pangenie.FAILED 2>/dev/null || echo "upstream PanGenie output missing or empty"; } > pangenie_var_plot.FAILED
        exit 0
    fi
    set +e
    python ${params.SCRIPT}/pangenie_var_plot.py $vcf0
    rc=\$?
    set -e
    if [ \$rc -ne 0 ] || [ ! -s pangenie_var_plot.png ]; then
        echo "pangenie_var_plot exited with status \$rc" > pangenie_var_plot.FAILED
        rm -f pangenie_var_plot.png
    fi
    exit 0
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
    tuple val(id), path("{pangenie_genotyping_biallelic.vcf.gz*,pangenie.FAILED}")

    tag "$id"
    publishDir "${params.outdir}/$id/", mode: 'link'

    stub:
    "touch pangenie_genotyping_biallelic.vcf.gz pangenie_genotyping_biallelic.vcf.gz.tbi"

    script:
    def bam_file = bam.first()
    def python = "/usr/local/app/miniconda3/bin/python"
    def py = "/usr/local/app/pangenie/pipelines/run-from-callset/scripts/convert-to-biallelic.py"
    """
    # Never fail the task: the final report must still be generated.
    # On failure write pangenie.FAILED + empty placeholder VCF; downstream plots then emit *.FAILED.
    set +e
    (
    set -e
        ${params.BIN}samtools fastq -@ ${task.cpus} -0 /dev/null -s /dev/null \
            -1 reads_r1.fq.gz -2 reads_r2.fq.gz $bam_file

        cat reads_r1.fq.gz reads_r2.fq.gz | gunzip > merge.fq

        PanGenie -f ${params.DB}/pangenie/HPRC_index -i merge.fq -o pangenie -j ${task.cpus} -t ${task.cpus}

        cat pangenie_genotyping.vcf | $python $py ${params.DB}/pangenie/cactus_filtered_ids_biallelic.vcf.gz | bgzip > pangenie_genotyping_biallelic.vcf.gz
        tabix pangenie_genotyping_biallelic.vcf.gz
    )
    rc=\$?
    set -e
    rm -f merge.fq pangenie_genotyping.vcf reads_r1.fq.gz reads_r2.fq.gz
    if [ \$rc -ne 0 ]; then
        echo "PanGenie exited with status \$rc" | tee pangenie.FAILED >&2
        : > pangenie_genotyping_biallelic.vcf.gz
        : > pangenie_genotyping_biallelic.vcf.gz.tbi
    fi
    exit 0
    """
}
