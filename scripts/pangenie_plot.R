require(RIdeogram)
args = (commandArgs(TRUE))

data(human_karyotype, package="RIdeogram")

contigs <- read.table("contigs.tsv", sep = "\t", header = T, stringsAsFactors = F)

SV <- read.table(args[1], sep = "\t", header = T, stringsAsFactors = F)

ideogram(karyotype = human_karyotype, overlaid = contigs, label = SV, label_type = "marker", colorset1 = c("#627A9D", "#77AAE4"))

convertSVG("chromosome.svg", device = "png")
