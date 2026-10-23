process specimmune {
    cpus params.cpu3
    memory params.specimmune_memory
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    when:
    !params.skip_specimmune && !params.demo && (params.ref == 'hg38' || params.ref.contains('GRCh38'))

    input:
    tuple val(id), path(bam)

    output:
    path 'specimmune_out'

    tag "$id"
    publishDir "${params.outdir}/report/$id/", mode: 'copy'

    script:
    def bams = bam instanceof List ? bam : [bam]
    bam = bams.first()
    """
    mkdir -p specimmune_out/HLA

    # Never fail the task: the final report must still be generated.
    # On failure write specimmune_out/HLA/FAILED (read by scripts/my_html.py) and exit 0.
    set +e
    (
    set -e
    if [ ! -x "${params.specimmune_env}/bin/python" ]; then
        echo "Missing SpecImmune Conda environment: --specimmune_env" >&2
        exit 2
    fi
    if [ ! -f "${params.specimmune_home}/scripts/main.py" ]; then
        echo "Missing SpecImmune source tree: --specimmune_home" >&2
        exit 2
    fi

    export PATH="${params.specimmune_env}/bin:\$PATH"

    # Feed SpecImmune only HLA-relevant reads (MHC region, chr6 MHC alt contigs, unmapped).
    # Whole-genome input made SpecImmune run >1 day and crash in visualization.py (>1M SVG elements).
    bai=""
    for f in ${bams.join(' ')}; do case "\$f" in *.bai|*.csi|*.crai) bai="\$f";; esac; done
    [ -n "\$bai" ] || ${params.BIN}samtools index -@ ${task.cpus} $bam
    regions="${params.specimmune_region}"
    for c in \$(${params.BIN}samtools idxstats $bam | cut -f1 | grep -E '${params.specimmune_alt_regex}' || true); do
        regions="\$regions \$c"
    done
    ${params.BIN}samtools view -u -F 0x900 -@ ${task.cpus} $bam \$regions '*' \\
        | ${params.BIN}samtools fastq -@ ${task.cpus} - | gzip -c > ${id}.specimmune.fastq.gz
    echo "SpecImmune input reads: \$(zcat ${id}.specimmune.fastq.gz | awk 'END{print NR/4}') (regions: \$regions)" >&2

    "${params.specimmune_env}/bin/python" "${params.specimmune_home}/scripts/main.py" \\
        -r ${id}.specimmune.fastq.gz -j ${task.cpus} -i HLA -n $id \\
        -o specimmune_out/HLA --db ${params.specimmune_db} \\
        --align_method_1 ${params.specimmune_align_method} -y ${params.specimmune_read_type}
    )
    rc=\$?
    set -e
    if [ \$rc -ne 0 ]; then
        echo "SpecImmune exited with status \$rc" | tee specimmune_out/HLA/FAILED >&2
    elif ! find specimmune_out/HLA -type f \\( -name '*.hap.alleles.txt' -o -name '*.final.type.result.formatted.txt' \\) -size +0 | grep -q .; then
        # SpecImmune sub-scripts can fail without a non-zero exit status (e.g. visualization.py traceback)
        echo "SpecImmune exited 0 but produced no HLA result file; see .command.err (last lines: \$(tail -n 3 .command.err 2>/dev/null | tr '\\n' ' '))" | tee specimmune_out/HLA/FAILED >&2
    fi
    exit 0
    """

    stub:
    "mkdir -p specimmune_out/HLA"
}
