process vep {
    
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    when: params.ref == 'hg38' || params.ref.contains('GRCh38')

    input:
    tuple val(id), path(phasedvcf)

    output:
    tuple val(id), path("${id}.vep.vcf.gz*"), emit: vcf
    tuple val(id), path("${id}.vep.vcf_summary.html"), emit: html

    tag "$id"
    publishDir "${params.outdir}/$id/annot/", pattern: "*vcf.gz*", mode: 'link'
    publishDir "${params.outdir}/report/$id/", pattern: "*html", mode: 'link'

    script:
    def vcf = phasedvcf.first()
    def ref = "${params.DB}/hg38/reference/hg38.fa"
    def db = "${params.DB}/hg38/"
    """
    set +u
    source /usr/local/app/miniconda3/bin/activate /usr/local/app/miniconda3/envs/vep

    vep -i $vcf -o ${id}.vep.vcf --cache --offline --dir_cache $db --fasta $ref --species homo_sapiens --assembly GRCh38 --cache_version 113 --vcf --fields Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,SYMBOL,dbNSF,AlphaMissense,UTRAnnotator,BIOTYPE,DISTANCE,FLAGS,VARIANT_CLASS,CLIN_SIG,AF

    bgzip -c ${id}.vep.vcf > ${id}.vep.vcf.gz
    tabix -p vcf ${id}.vep.vcf.gz
    """
}

process vep_data {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    
    input:
    tuple val(id), path(html)

    output:
    tuple val(id), path("*.{csv,png}")

    tag "$id"
    publishDir "${params.outdir}/report/$id/", mode: 'link'

    // cache false
    script:
    """
    python ${params.SCRIPT}/mk_vep_pic.py ${html}
    """
}