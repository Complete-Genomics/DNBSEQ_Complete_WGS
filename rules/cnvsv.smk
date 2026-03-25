# CNV and SV detection rules

# CNV detection using LFR-cnv
rule cnv_detection:
    input:
        bam = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam",
        vcf = config["outdir"] + "/{sample}/vcf/{sample}.{aligner}.{varcaller}.vcf.gz",
        hapblock = config["outdir"] + "/{sample}/phase/{sample}.{aligner}.{varcaller}.hapblock"
    output:
        cnv = config["outdir"] + "/{sample}/cnv/{sample}.CNV.result.xls"
    params:
        sample = "{sample}",
        aligner = "{aligner}",
        varcaller = "{varcaller}",
        ref = config.get("reference_name", "hg38"),
        script_dir = config["script_dir"]
    threads: config.get("threads", 8)
    container:
        config["container"]
    shell:
        """
        mkdir -p $(dirname {output.cnv})
        cd $(dirname {output.cnv})
        
        # Set reference-specific parameters
        if [ "{params.ref}" = "hs37d5" ]; then
            chr_param="-chr N"
            ref_param=""
        elif [ "{params.ref}" = "hg38" ]; then
            chr_param=""
            ref_param="-ref GRCH38"
        else
            chr_param=""
            ref_param=""
        fi
        
        # Set phase name pattern based on variant caller
        if [[ "{params.varcaller}" == *"dv"* ]]; then
            pname="{params.sample}.{params.aligner}.dv.XXX.hapblock"
        else
            pname="{params.sample}.{params.aligner}.gatk.XXX.hapblock"
        fi
        
        # Run LFR-cnv
        LFR-cnv -ncpu {threads} $ref_param \\
            -bam {input.bam} \\
            -vcf {input.vcf} \\
            -phase $(dirname {input.hapblock}) \\
            -pname $pname \\
            -sp 0.001 \\
            -out tmp \\
            -lcnv 1000 \\
            $chr_param
        
        mv tmp/ALL.200.format.cnv.1000.highconfidence {output.cnv}
        rm -rf tmp
        """

# SV detection using smoove
rule sv_detection:
    input:
        bam = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam",
        bai = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam.bai"
    output:
        vcf = config["outdir"] + "/{sample}/sv/{sample}.smoove.vcf.gz"
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
        cd $(dirname {output.vcf})
        
        # Run smoove for SV detection
        smoove call --outdir . --name {params.sample} --fasta {params.ref} -p {threads} {input.bam}
        
        # Rename output file
        mv {params.sample}-smoove.genotyped.vcf.gz {output.vcf}
        """
