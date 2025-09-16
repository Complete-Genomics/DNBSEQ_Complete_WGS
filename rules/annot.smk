# Annotation rules

# VEP annotation
rule vep_annotate:
    input:
        vcf = config["outdir"] + "/{sample}/phase/{sample}.{aligner}.{varcaller}.phased.vcf.gz"
    output:
        vep_vcf = config["outdir"] + "/{sample}/annot/{sample}.{aligner}.{varcaller}.vep.vcf.gz",
        vep_tbi = config["outdir"] + "/{sample}/annot/{sample}.{aligner}.{varcaller}.vep.vcf.gz.tbi"
    params:
        sample = "{sample}",
        aligner = "{aligner}",
        varcaller = "{varcaller}",
        ref = config["reference"],
        vep_cache = config["vep_cache_dir"]
    threads: config.get("threads", 4)
    container:
        config.get("vep_container", config["container"])
    shell:
        """
        mkdir -p $(dirname {output.vep_vcf})
        
        vep -i {input.vcf} \\
            -o temp_vep.vcf \\
            --cache --offline \\
            --dir_cache {params.vep_cache} \\
            --fasta {params.ref} \\
            --species homo_sapiens \\
            --assembly GRCh38 \\
            --cache_version 113 \\
            --vcf \\
            --fields Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,SYMBOL,dbNSF,AlphaMissense,UTRAnnotator,BIOTYPE,DISTANCE,FLAGS,VARIANT_CLASS,CLIN_SIG,AF
        
        bgzip -c temp_vep.vcf > {output.vep_vcf}
        tabix -p vcf {output.vep_vcf}
        rm temp_vep.vcf
        """
