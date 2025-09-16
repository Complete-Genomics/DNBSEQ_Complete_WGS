# HLA typing rules

# HLA typing with HLA-LA
rule hla_typing:
    input:
        bam = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam",
        bai = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam.bai"
    output:
        hla_out = directory(config["outdir"] + "/{sample}/hla/hlala_out")
    params:
        sample = "{sample}",
        aligner = "{aligner}",
        hla_db = config.get("hlala_db", "/path/to/hlala/PRG_MHC_GRCh38_withIMGT")
    threads: config.get("threads", 8)
    container:
        config["container"]
    shell:
        """
        mkdir -p $(dirname {output.hla_out})
        
        # Check if HLA reference is available
        if [ ! -e {params.hla_db}/knownReferences/hg38_corrected.txt ]; then
            samtools idxstats {input.bam} > tmp
            (
                echo -e "contigID\\tcontigLength\\tExtractCompleteContig\\tPartialExtraction_Start\\tPartialExtraction_Stop"
                head -n -1 tmp | awk -v OFS='\\t' '
                    {{
                        if ($1 == "chr6") {{
                            print $1, $2, 0, 28510120, 33480577
                        }} else {{
                            print $1, $2, 0
                        }}
                    }}'
                
                tail -n 1 tmp | awk -v OFS='\\t' '{{print $1, $2, 1}}'
            ) > hg38_corrected.txt
            mv hg38_corrected.txt {params.hla_db}/knownReferences/
        fi
        
        # Run HLA typing
        mkdir -p {output.hla_out}
        HLA-LA.pl --BAM {input.bam} --graph {params.hla_db} --sampleID {params.sample} --workingDir {output.hla_out} --maxThreads {threads}
        """
