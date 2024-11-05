suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(cowplot))

print(paste("Threads number is setting as", getDTthreads()))
option_list <- list(
  make_option(c("-t", "--table"), default = "NA12878_WES_v2/NA12878_WGS_v2_phased_possorted.barcode.stat.txt.gz",
              help = "long table of barcode - chromosome", dest = "filepath")
)
args <- parse_args(OptionParser(option_list = option_list))
filepath <- args$filepath

print("reading table file ...")
tab <- fread(filepath, sep = " ", sep2 = "\t", header = F)
tab <- data.frame(
  "num" = tab$V1,
  "chr" = gsub("(.+)\t(.+)", "\\1", tab$V2),
  "bcd" = gsub("(.+)\t(.+)", "\\2", tab$V2)
)
tab$chr <- gsub("(.+):(.+)", "\\1", tab$chr)

print("rearrange the table ...")
tab <- subset(tab, chr %in% paste0("chr", c(1:22, "X", "Y")))
print(paste("There is", nrow(tab), "records"))
tab <- recast(tab, bcd~chr)
rownames(tab) <- tab$bcd; tab$bcd <- NULL; tab[is.na(tab)] <- 0
tab <- as(as.matrix(tab), "dgCMatrix")

print("calculating common metric ...")
## fragment number
plot.data <- data.frame(
  row.names = rownames(tab),
  "readnum" = rowSums(tab),
  "shannon" = diversity(tab, index = "shannon"),
  "simpson" = diversity(tab, index = "simpson"),
  "mainpro" = apply(tab, 1, max)/rowSums(tab)
)

plot.data$readnumBin <- cut(plot.data$readnum,
                            breaks = c(0, 100, 2000, max(plot.data$readnum)))
plot.data$simpsonBin <- cut(plot.data$simpson,
                            breaks = c(0, 0.25, 0.6, 1), include.lowest = T)
plot.data$mainproBin <- cut(plot.data$simpson,
                            breaks = c(0, 0.45, 0.8, 1), include.lowest = T)

saveRDS(plot.data, file = gsub("txt.gz", "rds", filepath))

print("outputing Figures ...")
p1 <- ggplot(plot.data, aes(x = readnum)) +
  geom_density(fill = "steelblue", color = "steelblue") +
  xlim(0, 4000) +
  geom_vline(xintercept = c(100, 2000), color = "red") +
  theme_classic()

p2 <- ggplot(plot.data, aes(x = simpson)) +
  geom_density(fill = "steelblue", color = "steelblue") +
  geom_vline(xintercept = c(0.25, 0.6), color = "red") +
  theme_classic()

pdf(file = gsub("txt.gz", "density.pdf", filepath), width = 10, height = 5)
plot_grid(p1, p2)
dev.off()

