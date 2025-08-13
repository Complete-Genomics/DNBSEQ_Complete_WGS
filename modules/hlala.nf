process hlala {
    cpus params.cpu3
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    input:
    tuple val(id), path(bam)

    output:
    path "hlala_out"

    tag "$id"
    publishDir "${params.outdir}/report/$id/"

    script:
    bam = bam.first()
    """
    if [ ! -e ${params.DB}/hlala/PRG_MHC_GRCh38_withIMGT/knownReferences/hg38_corrected.txt ];then
        samtools idxstats $bam > tmp
        (
            echo -e "contigID\\tcontigLength\\tExtractCompleteContig\\tPartialExtraction_Start\\tPartialExtraction_Stop"
            head -n -1 tmp | awk -v OFS='\\t' '
                {
                    if (\$1 == "chr6") {
                        print \$1, \$2, 0, 28510120, 33480577
                    } else {
                        print \$1, \$2, 0
                    }
                }'
            
            tail -n 1 tmp | awk -v OFS='\\t' '{print \$1, \$2, 1}'
            ) > hg38_corrected.txt
        mv hg38_corrected.txt ${params.DB}/hlala/PRG_MHC_GRCh38_withIMGT/knownReferences/
    fi


    mkdir hlala_out
    HLA-LA.pl --BAM $bam --graph ${params.DB}/hlala/PRG_MHC_GRCh38_withIMGT --sampleID $id --workingDir hlala_out --maxThreads ${task.cpus}
    """
}