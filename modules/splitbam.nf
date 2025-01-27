// process splitBam {
  
//     cpus params.CPU0
//     memory params.MEM1 + "g"
//   clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
  

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
  clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
  

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
  # ${params.BIN}samtools view -h -F 0x400 ${bam} ${chr} \\
  #  | awk -v OFS='\\t' '{if(\$1~/#/){split(\$1,a,"#"); if(a[2]!~/0_0_0/){sub(/BX:Z:.*/, "BX:Z:"a[2]); print \$0} }else{print}}' - \\
  #  | ${params.BIN}samtools view -bhS - > ${id}.${aligner}.${chr}.bam

  # ${params.BIN}samtools view -h -F 0x400 $bam $chr \\
  # | perl -ne '\$p1=index(\$_,"#");\$p2=index(\$_,"\\t",\$p1);if (\$p1>0 && \$p2>0) {\$p3=rindex(\$_,"\\tBX:Z:");if (\$p3>0){substr(\$_,\$p3)="\\tBX:Z:".substr(\$_,\$p1+1,\$p2-\$p1-1)."\\n"}}print' | ${params.BIN}samtools view -bhS - > ${id}.${aligner}.${chr}.bam
  
  
  ${params.BIN}samtools view -h -F 0x400 $bam $chr \\
  | perl -ne '\$p1=index(\$_,"#");\$p2=index(\$_,"\\t",\$p1);if (\$p1>=0 && \$p2>\$p1) {\$bx_tag="BX:Z:".substr(\$_,\$p1+1,\$p2-\$p1-1); \$p3=rindex(\$_,"\tBX:Z:"); if (\$p3>=0){substr(\$_,\$p3)="\\t\$bx_tag\\n"} else {chomp; \$_.="\\t\$bx_tag\\n"}} print' | ${params.BIN}samtools view -bhS - > ${id}.${aligner}.${chr}.bam

  ${params.BIN}samtools index ${id}.${aligner}.${chr}.bam 
  """
  stub:
  "touch ${id}.${aligner}.${chr}.bam "
}

