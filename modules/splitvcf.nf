process split_vcf {
  
  cpus params.CPU0
  memory params.MEM0 + "g"
  clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
    
  input:
  val(aligner)
  val(varcaller)
  tuple val(id), path(vcf)
  each chr

  output:
  tuple val(id), val(chr), path("*.vcf.gz"), emit: eachvcf
  tuple val(id), path("*.vcf.gz"), emit: vcf

  tag "$id, $aligner, $varcaller"
  // publishDir "${params.outdir}/$id/align/alignsplit"

  script:
  def vcf = vcf.first()
  """
  zcat $vcf | awk -F'\\t' -v chr="$chr" '(\$0 ~ /^#/) || (\$1 == chr && \$10 ~ /^(1\\/1|0\\/1|1\\/2):/)' | ${params.BIN}bgzip > ${id}.${aligner}.${varcaller}.${chr}.vcf.gz
  """
}
