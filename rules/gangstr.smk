# GangSTR tandem repeat genotyping rules

# GangSTR genotyping
rule gangstr_genotyping:
    input:
        bam = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam",
        bai = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam.bai"
    output:
        vcf = config["outdir"] + "/{sample}/gangstr/GangSTR_out.vcf",
        insdata = config["outdir"] + "/{sample}/gangstr/GangSTR_out.insdata.tab"
    params:
        sample = "{sample}",
        aligner = "{aligner}",
        ref = config["reference"],
        gangstr_bed = config.get("gangstr_bed", "/path/to/GANGstr/hg38_ver17.bed")
    threads: config.get("threads", 8)
    container:
        config["container"]
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        
        # Run GangSTR for tandem repeat genotyping
        if [ -f "{params.gangstr_bed}" ]; then
            GangSTR --bam {input.bam} \\
                --ref {params.ref} \\
                --regions {params.gangstr_bed} \\
                --out $(dirname {output.vcf})/GangSTR_out \\
                --readlength {config[read_len]} \\
                --coverage 30 \\
                --insertmean 500 \\
                --insertsdev 50
        else
            echo "GangSTR BED file not found: {params.gangstr_bed}" > {output.vcf}
            touch {output.insdata}
        fi
        """
