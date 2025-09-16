# FASTQ processing rules

# Barcode splitting for stLFR data
rule barcode_split:
    input:
        fq1 = lambda wildcards: SAMPLES[wildcards.sample][wildcards.lane]['fqs'][0],
        fq2 = lambda wildcards: SAMPLES[wildcards.sample][wildcards.lane]['fqs'][1]
    output:
        split_fq1 = config["outdir"] + "/{sample}/split/{sample}_{lane}_1.fq.gz",
        split_fq2 = config["outdir"] + "/{sample}/split/{sample}_{lane}_2.fq.gz",
        split_log = config["outdir"] + "/{sample}/split/{sample}_{lane}.splitStat.xls"
    params:
        sample = "{sample}",
        lane = "{lane}",
        outdir = lambda wildcards: config["outdir"] + f"/{wildcards.sample}/split"
    threads: config.get("threads", 8)
    container:
        config["container"]
    shell:
        """
        mkdir -p {params.outdir}
        python {config[script_dir]}/stlfr2lariat_v3.py \\
            -1 {input.fq1} \\
            -2 {input.fq2} \\
            -o {params.outdir}/{params.sample}_{params.lane} \\
            -t {threads}
        """

# FASTQ quality control
rule fq_check:
    input:
        fq = lambda wildcards: SAMPLES[wildcards.sample][wildcards.lane]['fqs']
    output:
        report = config["outdir"] + "/{sample}/qc/{sample}_{lane}.fqcheck"
    params:
        sample = "{sample}",
        lane = "{lane}"
    container:
        config["container"]
    shell:
        """
        {config[script_dir]}/fqcheck/fqcheck_distribute.pl \\
            -f {input.fq[0]} \\
            -r {input.fq[1]} \\
            -o {output.report}
        """
