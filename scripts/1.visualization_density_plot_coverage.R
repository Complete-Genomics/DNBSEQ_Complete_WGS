# density plot of coverage

library(ggplot2)

args = (commandArgs(TRUE))
region_path = args[1]
density_plot = args[2]

df = data.frame()
txt = read.table(region_path)$V1
for(i in 1:length(txt)){

  file = txt[i]
  xdf = read.table(file)
  xdf$Data = strsplit(basename(file), "[.]")[[1]][1]
  df = rbind(df, xdf)
}
sample.col <- setNames(object = c("#2CA02CFF", "#1F77B4FF", "#1B191999", "#D62728FF", "#F79D1E99", "#FF7F0EFF"), nm = unique(df$Data))
pdf(density_plot, 12, 6)                          
p <- ggplot(df,
            aes(x = V5, color = Data)) +
            geom_density(bw = 1, linewidth = 1.2) +
            xlim(0, 60) +
            scale_color_manual(values = sample.col) +
            ylab("Density") + 
            xlab("Target Gene Coverage") +
            theme_classic() +
            theme(legend.position = "right",plot.title = element_text(hjust = 0.5)) +
            theme(axis.title.x =element_text(face = "bold", size=14), axis.title.y=element_text(face = "bold", size=14)) +
            theme(axis.text.x =element_text(size=12), axis.text.y=element_text(size=12)) + 
            theme(legend.title = element_text(face = "bold", size=14), legend.text = element_text(face = "bold", size=12))
print(p)                  
dev.off()



