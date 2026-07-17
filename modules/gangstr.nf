process gangstr {
    cpus params.CPU1
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    when: !params.demo && (params.ref == 'hg38' || params.ref.contains('GRCh38'))

    input:
    tuple val(id), path(bam)

    output:
    path "GangSTR_out*"

    tag "$id"
    publishDir "${params.outdir}/report/$id/", mode: 'copy'

    script:
    ref = params.ref.startsWith('/') ? params.ref : "${params.DB}/${params.ref}/reference/${params.ref}.fa"
    bed = "${params.DB}/GANGstr/hg38_ver17.bed"
    bam = bam.first()
    """
    GangSTR --bam $bam --ref $ref --regions $bed --out GangSTR_out
    """
    stub:
    "touch GangSTR_out.vcf GangSTR_out.samplestats.tab GangSTR_out.insdata.tab"
}