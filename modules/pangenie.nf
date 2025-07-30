process pangenie {
    cpus params.cpu3
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(r1), path(r2)

    output:
    path "pangenie_genotyping_biallelic.vcf"

    tag "$id"
    publishDir "${params.outdir}/report/$id/"

    script:
    def python = "/usr/local/app/miniconda3/bin/python"
    def py = "/usr/local/app/pangenie/pipelines/run-from-callset/scripts/convert-to-biallelic.py"
    """
    cat $r1 $r2 | gunzip > merge.fq
    PanGenie -f ${params.DB}/pangenie/HPRC_index -i merge.fq -o pangenie -j ${task.cpus} -t ${task.cpus}

    cat pangenie_genotyping.vcf | $python $py ${params.DB}/pangenie/cactus_filtered_ids_biallelic.vcf.gz > pangenie_genotyping_biallelic.vcf

    rm merge.fq
    """
}