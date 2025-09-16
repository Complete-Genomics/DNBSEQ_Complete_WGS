# Pangenie genotyping rules

# Pangenie genotyping
rule pangenie_genotyping:
    input:
        fq1 = config["outdir"] + "/{sample}/split/{sample}_{lane}_1.fq.gz",
        fq2 = config["outdir"] + "/{sample}/split/{sample}_{lane}_2.fq.gz"
    output:
        vcf = config["outdir"] + "/{sample}/pangenie/{sample}.pangenie.vcf"
    params:
        sample = "{sample}",
        pangenie_index = config.get("pangenie_index", "/path/to/pangenie/HPRC_index"),
        pangenie_biallelic = config.get("pangenie_biallelic", "/path/to/pangenie/cactus_filtered_ids_biallelic.vcf.gz")
    threads: config.get("threads", 16)
    container:
        config["container"]
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        
        # Run Pangenie genotyping
        if [ -f "{params.pangenie_index}" ] && [ -f "{params.pangenie_biallelic}" ]; then
            PanGenie -i {params.pangenie_index} \\
                -r {input.fq1},{input.fq2} \\
                -v {params.pangenie_biallelic} \\
                -o {params.sample} \\
                -t {threads}
            
            mv {params.sample}_genotyping.vcf {output.vcf}
        else
            echo "Pangenie index or biallelic VCF not found" > {output.vcf}
        fi
        """
