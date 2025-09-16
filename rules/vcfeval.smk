# VCF evaluation rules

# VCF statistics
rule vcf_stats:
    input:
        vcf = config["outdir"] + "/{sample}/vcf/{sample}.{aligner}.{varcaller}.vcf.gz"
    output:
        stats = config["outdir"] + "/{sample}/vcfeval/{sample}.{aligner}.{varcaller}.stats"
    params:
        sample = "{sample}",
        aligner = "{aligner}",
        varcaller = "{varcaller}"
    threads: config.get("threads", 4)
    container:
        config["container"]
    shell:
        """
        mkdir -p $(dirname {output.stats})
        
        # Calculate VCF statistics
        bcftools stats {input.vcf} > {output.stats}
        """

# VCF evaluation against benchmark
rule vcf_benchmark:
    input:
        vcf = config["outdir"] + "/{sample}/vcf/{sample}.{aligner}.{varcaller}.vcf.gz",
        benchmark_vcf = config.get("benchmark_vcf", "")
    output:
        eval_report = config["outdir"] + "/{sample}/vcfeval/{sample}.{aligner}.{varcaller}.eval.txt"
    params:
        sample = "{sample}",
        aligner = "{aligner}",
        varcaller = "{varcaller}"
    threads: config.get("threads", 4)
    container:
        config["container"]
    shell:
        """
        mkdir -p $(dirname {output.eval_report})
        
        if [ -n "{input.benchmark_vcf}" ] && [ -f "{input.benchmark_vcf}" ]; then
            # Run VCF comparison
            bcftools isec -p vcf_comparison {input.vcf} {input.benchmark_vcf}
            echo "VCF evaluation completed" > {output.eval_report}
        else
            echo "No benchmark VCF provided" > {output.eval_report}
        fi
        """
