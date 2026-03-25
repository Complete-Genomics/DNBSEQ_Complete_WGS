# Alignment rules

# Alignment with BWA/VG
rule align:
    input:
        fq1 = config["outdir"] + "/{sample}/split/{sample}_{lane}_1.fq.gz",
        fq2 = config["outdir"] + "/{sample}/split/{sample}_{lane}_2.fq.gz"
    output:
        bam = config["outdir"] + "/{sample}/align/{sample}_{lane}.{aligner}.bam",
        bai = config["outdir"] + "/{sample}/align/{sample}_{lane}.{aligner}.bam.bai"
    params:
        sample = "{sample}",
        lane = "{lane}",
        aligner = "{aligner}",
        ref = config["reference"]
    threads: config.get("threads", 16)
    container:
        config["container"]
    shell:
        """
        if [ "{params.aligner}" = "bwa" ]; then
            bwa mem -t {threads} -M \\
                -R '@RG\\tID:{params.sample}_{params.lane}\\tSM:{params.sample}\\tLB:{params.sample}\\tPL:DNBSEQ' \\
                {params.ref} {input.fq1} {input.fq2} | \\
            samtools sort -@ {threads} -o {output.bam} -
        elif [ "{params.aligner}" = "vg" ]; then
            vg giraffe -t {threads} \\
                -R '@RG\\tID:{params.sample}_{params.lane}\\tSM:{params.sample}\\tLB:{params.sample}\\tPL:DNBSEQ' \\
                -x {params.ref}.xg -g {params.ref}.gcsa \\
                -f {input.fq1} -f {input.fq2} | \\
            samtools sort -@ {threads} -o {output.bam} -
        fi
        samtools index {output.bam}
        """

# Merge BAMs per sample
rule merge_bams:
    input:
        bams = lambda wildcards: expand(
            config["outdir"] + "/{sample}/align/{sample}_{lane}.{aligner}.bam",
            sample=wildcards.sample,
            lane=SAMPLES[wildcards.sample].keys(),
            aligner=wildcards.aligner
        )
    output:
        merged_bam = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam",
        merged_bai = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam.bai"
    params:
        sample = "{sample}",
        aligner = "{aligner}"
    threads: config.get("threads", 8)
    container:
        config["container"]
    shell:
        """
        if [ $(echo {input.bams} | wc -w) -eq 1 ]; then
            ln -sf $(realpath {input.bams}) {output.merged_bam}
            ln -sf $(realpath {input.bams}).bai {output.merged_bai}
        else
            samtools merge -@ {threads} {output.merged_bam} {input.bams}
            samtools index {output.merged_bam}
        fi
        """
