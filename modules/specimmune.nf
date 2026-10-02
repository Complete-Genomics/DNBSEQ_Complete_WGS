process specimmune {
    cpus params.cpu3
    memory params.specimmune_memory
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    when:
    !params.skip_specimmune && !params.demo && (params.ref == 'hg38' || params.ref.contains('GRCh38'))

    input:
    tuple val(id), path(bam)

    output:
    path 'specimmune_out'

    tag "$id"
    publishDir "${params.outdir}/report/$id/", mode: 'copy'

    script:
    bam = bam instanceof List ? bam.first() : bam
    """
    mkdir -p specimmune_out/HLA
    ${params.BIN}samtools fastq -@ ${task.cpus} -F 0x900 $bam | gzip -c > ${id}.specimmune.fastq.gz

    python3 /SpecImmune/scripts/main.py \\
        -r ${id}.specimmune.fastq.gz -j ${task.cpus} -i HLA -n $id \\
        -o specimmune_out/HLA --db ${params.specimmune_db} \\
        --align_method_1 ${params.specimmune_align_method} -y ${params.specimmune_read_type}
    """

    stub:
    "mkdir -p specimmune_out/HLA"
}
