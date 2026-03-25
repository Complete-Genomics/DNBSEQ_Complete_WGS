# BAM filtering and processing rules

# BAM quality filtering
rule filter_bam:
    input:
        bam = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam",
        bai = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam.bai"
    output:
        filtered_bam = config["outdir"] + "/{sample}/align/{sample}.{aligner}.filtered.bam",
        filtered_bai = config["outdir"] + "/{sample}/align/{sample}.{aligner}.filtered.bam.bai"
    params:
        sample = "{sample}",
        aligner = "{aligner}",
        mapq_threshold = config.get("mapq_threshold", 20)
    threads: config.get("threads", 8)
    container:
        config["container"]
    shell:
        """
        samtools view -@ {threads} -b -q {params.mapq_threshold} -F 1804 {input.bam} > {output.filtered_bam}
        samtools index {output.filtered_bam}
        """

# Mark duplicates
rule mark_duplicates:
    input:
        bam = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam"
    output:
        marked_bam = config["outdir"] + "/{sample}/align/{sample}.{aligner}.marked.bam",
        marked_bai = config["outdir"] + "/{sample}/align/{sample}.{aligner}.marked.bam.bai",
        metrics = config["outdir"] + "/{sample}/align/{sample}.{aligner}.markdup.metrics"
    params:
        sample = "{sample}",
        aligner = "{aligner}"
    threads: config.get("threads", 8)
    container:
        config["container"]
    shell:
        """
        picard MarkDuplicates \\
            INPUT={input.bam} \\
            OUTPUT={output.marked_bam} \\
            METRICS_FILE={output.metrics} \\
            CREATE_INDEX=true \\
            VALIDATION_STRINGENCY=LENIENT
        """
