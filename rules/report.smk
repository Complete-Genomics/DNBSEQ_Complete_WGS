# Report generation rules

# Generate sample report
rule sample_report:
    input:
        vcf = config["outdir"] + "/{sample}/annot/{sample}.{aligner}.{varcaller}.vep.vcf.gz",
        phase_report = config["outdir"] + "/{sample}/phase/{sample}.{aligner}.{varcaller}.phase.report",
        chromosome_png = config["outdir"] + "/{sample}/phase/chromosome.png"
    output:
        report = config["outdir"] + "/{sample}/report/{sample}.{aligner}.{varcaller}.report"
    params:
        sample = "{sample}",
        aligner = "{aligner}",
        varcaller = "{varcaller}",
        script_dir = config["script_dir"]
    container:
        config["container"]
    shell:
        """
        mkdir -p $(dirname {output.report})
        
        # Calculate variant statistics
        snp=$(bcftools view -v snps {input.vcf} | grep -v \\# | wc -l)
        indel=$(bcftools view -v indels {input.vcf} | grep -v \\# | wc -l)
        hetsnp=$(bcftools view -v snps -i 'GT="0/1" || GT="1|0" || GT="0|1"' {input.vcf} | grep -v \\# | wc -l)
        hetindel=$(bcftools view -v indels -i 'GT="0/1" || GT="1|0" || GT="0|1"' {input.vcf} | grep -v \\# | wc -l)
        phasedhetsnp=$(bcftools view -v snps -i 'GT="1|0" || GT="0|1"' {input.vcf} | grep -v \\# | wc -l)
        phasedhetindel=$(bcftools view -v indels -i 'GT="1|0" || GT="0|1"' {input.vcf} | grep -v \\# | wc -l)
        
        echo -e "$snp\\t$indel\\t$hetsnp\\t$hetindel\\t$phasedhetsnp\\t$phasedhetindel" > $(dirname {output.report})/varstat
        
        # Generate report using existing script
        python {params.script_dir}/report.py ref {params.sample} {input.vcf} \\
            dummy_lfr dummy_depth {input.phase_report} > {output.report}
        """

# Generate HTML report
rule html_report:
    input:
        reports = expand(
            config["outdir"] + "/{sample}/report/{sample}.{aligner}.{varcaller}.report",
            sample=SAMPLE_IDS,
            aligner=ALIGNERS,
            varcaller=VARCALLERS
        ),
        chromosome_pngs = expand(
            config["outdir"] + "/{sample}/phase/chromosome.png",
            sample=SAMPLE_IDS
        )
    output:
        html_reports = expand(
            config["outdir"] + "/{sample}/report/{sample}_report.html",
            sample=SAMPLE_IDS
        )
    params:
        outdir = config["outdir"] + "/report/",
        script_dir = config["script_dir"]
    container:
        config["container"]
    shell:
        """
        mkdir -p {params.outdir}
        python {params.script_dir}/my_html.py {params.outdir}
        """

# Combine final reports
rule combine_reports:
    input:
        reports = expand(
            config["outdir"] + "/{sample}/report/{sample}.{aligner}.{varcaller}.report",
            sample=SAMPLE_IDS,
            aligner=ALIGNERS,
            varcaller=VARCALLERS
        )
    output:
        csv = config["outdir"] + "/report.csv"
    shell:
        """
        paste {input.reports} | awk -F'\\t' '{{
            line = $1
            for (i = 2; i <= NF; i += 2) {{
                line = line "," $i
            }}
            print line
        }}' > {output.csv}
        """

# Generate FASTQ QC report
rule fqc_report:
    input:
        fq_reports = expand(
            config["outdir"] + "/{sample}/qc/{sample}_{lane}.fqcheck",
            sample=SAMPLE_IDS,
            lane=[lane for sample in SAMPLE_IDS for lane in SAMPLES[sample].keys()]
        )
    output:
        fqc_csv = config["outdir"] + "/report_fqc.csv"
    params:
        script_dir = config["script_dir"]
    container:
        config["container"]
    shell:
        """
        python {params.script_dir}/fqc.py > {output.fqc_csv}
        """
