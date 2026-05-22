// process splitBam {
  
//     cpus params.CPU0
//     memory params.MEM1 + "g"
//   clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
  

//   input:
//   tuple val(id), path(bam)
//   each chr

//   output:
//   tuple val(id), val(chr), path("${id}.lariat0.${chr}.bam*"), emit: eachbam

//   tag "$id split lariatbam for markdup"
  
//   publishDir "${params.outdir}/$id/align/alignsplit"

//   script:
//   def bam = bam.first()
//   """
//   ${params.BIN}samtools view -bh ${bam} ${chr} -o ${id}.lariat0.${chr}.bam
//   ${params.BIN}samtools index ${id}.lariat0.${chr}.bam 
//   """
// }

process splitBam4phasing {
  
  cpus params.CPU0
  memory params.MEM0 + "g"
  clusterOptions = params.clusterOptions.replace('CPUS', cpus.toString()).replace('MEMORY', memory.toString()).replace('QUEUE', params.queue)
  

  input:
  val(aligner)
  tuple val(id), path(bam)
  each chr

  output:
  tuple val(id), val(chr), path("${id}.${aligner}.${chr}.bam*"), emit: eachbam

  tag "$id,$aligner"
  
  // publishDir "${params.outdir}/$id/align/alignsplit"

  script:
  def bam = bam.first()
  def prefix = "${bam.getBaseName()}"
  """
  # Auto-detect library type by presence of BX:Z: tag in the first 100 reads.
  #   SE600: BX:Z: already present; strip BX:Z:0_0_0 (no valid barcode) before passing through.
  #   PE150 lariat: no BX tag; derive BX from #umi in QNAME via perl.
  n_bx=\$(${params.BIN}samtools view $bam $chr 2>/dev/null | head -100 | grep -c "BX:Z:" || true)
  if [ "\$n_bx" -gt 0 ]; then
      # SE600 path — strip BX:Z:0_0_0 (reads with no valid barcode) to avoid fake linked-fragments in extractHAIRS
      ${params.BIN}samtools view -h -F 0x400 $bam $chr \\
      | awk '{if(/^@/){print; next} gsub(/\\tBX:Z:0_0_0/,""); print}' \\
      | ${params.BIN}samtools view -bh -o ${id}.${aligner}.${chr}.bam
  else
      # PE150 lariat path — read name has #umi, perl writes BX from it
      ${params.BIN}samtools view -h -F 0x400 $bam $chr \\
      | perl -ne '\$p1=index(\$_,"#");\$p2=index(\$_,"\\t",\$p1);if (\$p1>=0 && \$p2>\$p1) {\$bx_tag="BX:Z:".substr(\$_,\$p1+1,\$p2-\$p1-1); \$p3=rindex(\$_,"\\tBX:Z:"); if (\$p3>=0){substr(\$_,\$p3)="\\t\$bx_tag\\n"} else {chomp; \$_.="\\t\$bx_tag\\n"}} print' \\
      | ${params.BIN}samtools view -bhS - > ${id}.${aligner}.${chr}.bam
  fi
  ${params.BIN}samtools index ${id}.${aligner}.${chr}.bam
  """
  stub:
  "touch ${id}.${aligner}.${chr}.bam "
}

