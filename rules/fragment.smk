# Fragment analysis rules

# Fragment length analysis
rule fragment_analysis:
    input:
        bam = config["outdir"] + "/{sample}/align/{sample}.{aligner}.merged.bam"
    output:
        fragment_stats = config["outdir"] + "/{sample}/fragment/{sample}.{aligner}.fragment.stats"
    params:
        sample = "{sample}",
        aligner = "{aligner}"
    threads: config.get("threads", 4)
    container:
        config["container"]
    shell:
        """
        mkdir -p $(dirname {output.fragment_stats})
        
        # Calculate fragment size statistics
        samtools view -f 2 {input.bam} | \\
        awk '{{if($9>0) print $9}}' | \\
        sort -n | \\
        awk '{{
            count++
            sum+=$1
            if(count==1) min=$1
            max=$1
        }}
        END{{
            print "Sample: {params.sample}"
            print "Aligner: {params.aligner}"
            print "Total fragments: " count
            print "Mean fragment size: " sum/count
            print "Min fragment size: " min
            print "Max fragment size: " max
        }}' > {output.fragment_stats}
        """
