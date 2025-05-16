require(RIdeogram)

args<-commandArgs(TRUE)
genome <- args[1]

human_karyotype <- read.table(genome, sep = "\t", header = T, stringsAsFactors = F)
contigs <- read.table("contigs.tsv", sep = "\t", header = T, stringsAsFactors = F)
ideogram(karyotype = human_karyotype, overlaid = contigs, colorset1 = c("#627A9D", "#77AAE4"), Lx = -200, Ly = -200)
convertSVG("chromosome.svg", device = "png",width=10, height=7)