# Quality control rules

# Basic quality control statistics
rule qc_stats:
    input:
        fq1 = config["outdir"] + "/{sample}/split/{sample}_{lane}_1.fq.gz",
        fq2 = config["outdir"] + "/{sample}/split/{sample}_{lane}_2.fq.gz"
    output:
        stats = config["outdir"] + "/{sample}/qc/{sample}_{lane}.qc.stats"
    params:
        sample = "{sample}",
        lane = "{lane}"
    threads: config.get("threads", 4)
    container:
        config["container"]
    shell:
        """
        mkdir -p $(dirname {output.stats})
        
        # Basic QC statistics
        echo "Sample: {params.sample}" > {output.stats}
        echo "Lane: {params.lane}" >> {output.stats}
        echo "Read1 file: {input.fq1}" >> {output.stats}
        echo "Read2 file: {input.fq2}" >> {output.stats}
        
        # Count reads
        r1_reads=$(zcat {input.fq1} | wc -l | awk '{{print $1/4}}')
        r2_reads=$(zcat {input.fq2} | wc -l | awk '{{print $1/4}}')
        echo "Read1 count: $r1_reads" >> {output.stats}
        echo "Read2 count: $r2_reads" >> {output.stats}
        """
