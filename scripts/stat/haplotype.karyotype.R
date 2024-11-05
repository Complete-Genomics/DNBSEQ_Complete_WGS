
args<-commandArgs(TRUE)

genome <- args[1]
band <- args[2]
out <- args[3]
lib <- args[4]

library(regioneR)
library(karyoploteR)

pdf(out)
custom.genome <- toGRanges(genome)
custom.cytobands <- toGRanges(band)

pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight = 50
pp$data1height = 20

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, plot.params = pp)

dev.off()
