process avi {
    cpus params.CPU0
    memory params.MEM1 + "g"
    clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)

    when: !params.skip_avi && (params.ref == 'hg38' || params.ref.contains('GRCh38'))

    input:
    tuple val(id), path(vcf)    // VEP-annotated phased VCF (+ .tbi), so the top table carries gene/consequence

    output:
    tuple val(id), path("${id}.avi.vcf.gz*"), optional: true, emit: vcf
    path "avi_*", emit: report  // avi_top.tsv + avi_summary.tsv, or avi_top.FAILED

    tag "$id"
    publishDir "${params.outdir}/$id/annot/", pattern: "*.avi.vcf.gz*", mode: 'link'
    publishDir "${params.outdir}/report/$id/", pattern: "avi_*", mode: 'copy'

    stub:
    "touch ${id}.avi.vcf.gz ${id}.avi.vcf.gz.tbi avi_top.tsv avi_summary.tsv"

    script:
    def vcf0 = vcf instanceof List ? vcf.first() : vcf
    // prepared once with scripts/avi_prep.sh: CHROM POS REF ALT PHRED, bgzip + tabix
    def table = params.avi_table ?: "${params.DB}/hg38/alphagenome/alphagenome_avi_phred.tsv.gz"
    """
    # Soft-fail: write avi_top.FAILED (read by scripts/my_html.py) instead of failing the task.
    set +e
    (
    set -e
    if [ ! -s "${table}" ] || [ ! -s "${table}.tbi" ]; then
        echo "AVI table or its .tbi not found: ${table} (build it with scripts/avi_prep.sh)" >&2
        exit 2
    fi
    echo '##INFO=<ID=AVI,Number=A,Type=Float,Description="AlphaGenome variant impact score (PHRED), SNVs only">' > avi.hdr
    bcftools annotate -a ${table} -h avi.hdr -c CHROM,POS,REF,ALT,AVI -Oz -o ${id}.avi.vcf.gz ${vcf0}
    tabix -f -p vcf ${id}.avi.vcf.gz
    bcftools view ${id}.avi.vcf.gz | ${params.BIN}python3 ${params.SCRIPT}/avi_top.py ${params.avi_top_n} avi_top.tsv avi_summary.tsv
    ) 2> avi.err
    rc=\$?
    set -e
    cat avi.err >&2
    if [ \$rc -ne 0 ] || [ ! -s avi_top.tsv ]; then
        { echo "AVI exited with status \$rc:"; tail -n 3 avi.err; } | tr '\\n' ' ' > avi_top.FAILED
        rm -f avi_top.tsv avi_summary.tsv
    fi
    exit 0
    """
}
