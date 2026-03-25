# Gene depth analysis rules

# Gene coverage analysis
rule gene_coverage:
    input:
        bam = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam",
        bai = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam.bai"
    output:
        coverage = config["outdir"] + "/{sample}/genedepth/{sample}.{aligner}.gene.coverage"
    params:
        sample = "{sample}",
        aligner = "{aligner}",
        gene_bed = config.get("gene_bed", "genes.bed")
    threads: config.get("threads", 4)
    container:
        config["container"]
    shell:
        """
        mkdir -p $(dirname {output.coverage})
        
        # Calculate gene coverage
        if [ -f "{params.gene_bed}" ]; then
            bedtools coverage -a {params.gene_bed} -b {input.bam} > {output.coverage}
        else
            echo "Gene BED file not found: {params.gene_bed}" > {output.coverage}
        fi
        """

# Depth statistics
rule depth_stats:
    input:
        bam = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam"
    output:
        depth = config["outdir"] + "/{sample}/genedepth/{sample}.{aligner}.depth.stats"
    params:
        sample = "{sample}",
        aligner = "{aligner}"
    threads: config.get("threads", 4)
    container:
        config["container"]
    shell:
        """
        mkdir -p $(dirname {output.depth})
        
        # Calculate depth statistics
        samtools depth {input.bam} | \\
        awk '{{
            count++
            sum+=$3
            if($3>0) covered++
        }}
        END{{
            print "Sample: {params.sample}"
            print "Aligner: {params.aligner}"
            print "Total positions: " count
            print "Covered positions: " covered
            print "Coverage percentage: " (covered/count)*100 "%"
            print "Mean depth: " sum/count
        }}' > {output.depth}
        """
