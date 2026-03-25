# Phasing rules

# Phasing with HapCUT2
rule hapcut2_phase:
    input:
        vcf = config["outdir"] + "/{sample}/vcf/{sample}.{aligner}.{varcaller}.vcf.gz",
        bam = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam"
    output:
        hapblock = config["outdir"] + "/{sample}/phase/{sample}.{aligner}.{varcaller}.hapblock",
        phased_vcf = config["outdir"] + "/{sample}/phase/{sample}.{aligner}.{varcaller}.phased.vcf.gz",
        phase_report = config["outdir"] + "/{sample}/phase/{sample}.{aligner}.{varcaller}.phase.report"
    params:
        sample = "{sample}",
        aligner = "{aligner}",
        varcaller = "{varcaller}"
    threads: config.get("threads", 8)
    container:
        config["container"]
    shell:
        """
        mkdir -p $(dirname {output.hapblock})
        
        # Extract haplotype-informative reads
        extractHAIRS --bam {input.bam} --VCF {input.vcf} --out temp_fragment_file
        
        # Phase variants
        HAPCUT2 --fragments temp_fragment_file --VCF {input.vcf} --output {output.hapblock}
        
        # Generate phased VCF
        python {config[script_dir]}/phase.py {input.vcf} {output.hapblock} > {output.phased_vcf}
        
        # Generate phase statistics
        python {config[script_dir]}/calculate_haplotype_statistics.py \\
            {output.hapblock} > {output.phase_report}
        
        rm -f temp_fragment_file
        """

# Generate chromosome ideogram
rule ideogram:
    input:
        hapblock = config["outdir"] + "/{sample}/phase/{sample}.{aligner}.{varcaller}.hapblock"
    output:
        png = config["outdir"] + "/{sample}/phase/chromosome.png",
        svg = config["outdir"] + "/{sample}/phase/chromosome.svg"
    params:
        sample = "{sample}",
        script_dir = config["script_dir"]
    container:
        config["container"]
    shell:
        """
        cd $(dirname {output.png})
        
        # Convert hapblock to ideogram format (using modified band.py with SPAN filter)
        python {params.script_dir}/band.py {input.hapblock}
        
        # Generate chromosome ideogram
        Rscript {params.script_dir}/ideogram.R {params.script_dir}/hg38.karyotype
        
        # Crop and convert to PNG
        convert -crop 1250x1250+900+200 chromosome.png tmp
        mv tmp chromosome.png
        """
