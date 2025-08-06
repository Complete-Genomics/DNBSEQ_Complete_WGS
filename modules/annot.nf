process vep {
    
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
    
    input:
    tuple val(id), path(phasedvcf)

    output:
    tuple val(id), path("${id}.vep.vcf.gz*")

    tag "$id"
    publishDir "${params.outdir}/$id/annot/", mode: 'link'

    script:
    def vcf = phasedvcf.first()
    def ref = "${params.DB}/hg38/reference/hg38.fa"
    def db = "${params.DB}/hg38/"
    """
    set +u
    source /usr/local/miniconda3/bin/activate /usr/local/miniconda3/envs/vep

    vep -i $vcf -o ${id}.vep.vcf --cache --offline --dir_cache $db --fasta $ref --species homo_sapiens --assembly GRCh38 --cache_version 113 --vcf --fields Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,SYMBOL,dbNSF,AlphaMissense,UTRAnnotator

    bgzip -c ${id}.vep.vcf > ${id}.vep.vcf.gz
    tabix -p vcf ${id}.vep.vcf.gz
    """
}