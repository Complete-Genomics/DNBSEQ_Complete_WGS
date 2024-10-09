process mapq {
  
  cpus params.CPU1
  memory params.MEM1 + "g"
  clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"
  
  input:
  tuple val(id), path(bam)

  output:
  tuple val(id), path("${id}.pf.mapqlt${params.pfmapq}.bam*")

  tag "$id"  
  publishDir "${params.outdir}/$id/align"

  script:
  def bam = bam.first()
  if (params.pfmapq != 0) {
  	"""
  	${params.BIN}samtools view -b -q ${params.pfmapq} $bam -o ${id}.pf.mapqlt${params.pfmapq}.bam
  	${params.BIN}samtools index ${id}.pf.mapqlt${params.pfmapq}.bam
  	"""
  } else {
	"""
	cp $bam ${id}.pf.mapqlt${params.pfmapq}.bam
	cp ${bam}.bai ${id}.pf.mapqlt${params.pfmapq}.bam.bai
	"""
  }
  stub:
  "touch ${id}.pf.mapqlt${params.pfmapq}.bam "
}

