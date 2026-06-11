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
  # Parse BX from QNAME '#umi' field for both PE150 and SE600.
  # Safety: search for '#' only within QNAME (before the first tab) to prevent
  # matching '#' characters that appear in QUAL scores, which would corrupt the line.
  # Reads with barcode 0_0_0 (no valid bead barcode) skip BX tagging so they are
  # treated as unlinked singles by extractHAIRS, not grouped into a fake fragment.
  ${params.BIN}samtools view -h -F 0x400 $bam $chr \\
  | perl -ne '
      if (/^@/) { print; next }
      \$p0 = index(\$_, "\\t");
      \$p1 = index(\$_, "#");
      if (\$p1 >= 0 && \$p1 < \$p0) {
          \$p2 = index(\$_, "\\t", \$p1);
          \$bc = substr(\$_, \$p1+1, \$p2-\$p1-1);
          \$bc =~ s/\\/[12]\$//;
          if (\$bc ne "0_0_0") {
              \$bx_tag = "BX:Z:\$bc";
              \$p3 = rindex(\$_, "\\tBX:Z:");
              if (\$p3 >= 0) { substr(\$_, \$p3) = "\\t\$bx_tag\\n" }
              else           { chomp; \$_ .= "\\t\$bx_tag\\n" }
          }
      }
      print' \\
  | ${params.BIN}samtools view -bhS - > ${id}.${aligner}.${chr}.bam
  ${params.BIN}samtools index ${id}.${aligner}.${chr}.bam
  """
  stub:
  "touch ${id}.${aligner}.${chr}.bam "
}

