# Variant calling rules

# Variant calling with DeepVariant
rule deepvariant:
    input:
        bam = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam",
        bai = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam.bai"
    output:
        vcf = config["outdir"] + "/{sample}/vcf/{sample}.{aligner}.{varcaller}.vcf.gz",
        tbi = config["outdir"] + "/{sample}/vcf/{sample}.{aligner}.{varcaller}.vcf.gz.tbi"
    params:
        sample = "{sample}",
        aligner = "{aligner}",
        varcaller = "{varcaller}",
        ref = config["reference"],
        model = config.get("dv_model", "WGS")
    threads: config.get("threads", 16)
    container:
        config.get("deepvariant_container", config["container"])
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        /opt/deepvariant/bin/run_deepvariant \\
            --model_type={params.model} \\
            --ref={params.ref} \\
            --reads={input.bam} \\
            --output_vcf={output.vcf} \\
            --num_shards={threads}
        """

# Variant calling with GATK HaplotypeCaller
rule gatk_haplotypecaller:
    input:
        bam = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam",
        bai = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam.bai"
    output:
        vcf = config["outdir"] + "/{sample}/vcf/{sample}.{aligner}.gatk.vcf.gz",
        tbi = config["outdir"] + "/{sample}/vcf/{sample}.{aligner}.gatk.vcf.gz.tbi"
    params:
        sample = "{sample}",
        aligner = "{aligner}",
        ref = config["reference"]
    threads: config.get("threads", 8)
    container:
        config["container"]
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        gatk HaplotypeCaller \\
            -R {params.ref} \\
            -I {input.bam} \\
            -O temp.vcf.gz \\
            --native-pair-hmm-threads {threads}
        
        mv temp.vcf.gz {output.vcf}
        mv temp.vcf.gz.tbi {output.tbi}
        """
