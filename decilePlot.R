#!/usr/bin/env Rscript

source("D:/PCa_Discovery/R_scripts/mqlib_prot.R")
theme_set(plot_theme())

#setwd("D:/PCa_Discovery/Data_analysis_files/CellLine/txt_190918_pca_wcl/")

# read in pg df
pg <- pgRead("R_tableOutput/lfq_df.txt")

# get intensities
int <- pg[,grep("^lfq", names(pg))]

int <- int[,-grep("Blank", names(int))]

# remove all 0 columns
int <- int[,colSums(is.na(int)) < nrow(int)]

# number of observations
num.obs <- apply(int, 1, function(x) sum(!is.na(x)))

median <- apply(int, 1, function(x) median(x, na.rm = T))

# split proteins into deciles by median intensity
decile <- ntile(-median, 10)

# bind df
df <- cbind.data.frame(pg$genes, num.obs, median, decile)

# allCounts

exo.markers <- c("CD9", "CD63", "CD81", "PDCD6IP", "TSG101", "FLOT1", "FLOT2", "ANXA2", "ARRDC1")
prostate.markers <- c("KLK3", "ACPP", "TGM4", "FOLH1", "KLK2", "MSMB", "NKX3.1", "NEFH")

markers <- c("UMOD", "ALB", "APOA1", exo.markers, prostate.markers)

label <- df[rownames(df) %in% markers,]

protcount <- ggplot(df, aes(num.obs))+
  geom_histogram(binwidth = 1, color = "black", fill = "white")+
  scale_x_reverse()+
  labs(x = NULL, y = "Protein Count")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20))
  
#svg(file = paste(dir2, "ProteinAbundanceCurve_bottomPlot.svg", sep = ""), height = 5, width = 7)
scurve <- ggplot(df, aes(num.obs, median, color = factor(decile)))+
  geom_point(shape = 1, alpha = 0.6, size = 1)+
  scale_y_log10()+
  scale_color_brewer(palette = "Spectral")+
  scale_x_reverse()+
  labs(x = "Number of observations", y = "Median Intensity", color = "Decile")+
  geom_point(data = label, aes(num.obs, median), color = ifelse(rownames(label) %in% exo.markers, "blue", "black"), shape = 16, size = 0.8)+
  geom_text_repel(data = label,
                  label = rownames(label),
                  nudge_x = ifelse(rownames(label) %in% exo.markers, 10, 20),
                  color = ifelse(rownames(label) %in% exo.markers, "blue", ifelse(rownames(label) %in% prostate.markers, "black", "red")),
                  size = 5,
                  segment.alpha = 0.5,
                  force = 3)+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))
#dev.off()


svg(file = "R_figures/ProteinAbundanceCurve_legend.svg")
ggarrange(scurve, protcount, nrow = 2, align = "hv", heights = c(1, 0.5), common.legend = TRUE)
dev.off()

rm(list = ls())
q()
