#!/usr/bin/env Rscript

source("D:/PCa_Discovery/R_scripts/mqlib_prot.R")
theme_set(plot_theme())

setwd("D:/PCa_Discovery/Data_analysis_files/CellLine/txt_190918_pca_exo/")

# read in pg df

pg <- pgRead("R_tableOutput/pep_lfq_df.txt")

# get intensities
int <- pg[,grep("^lfq", names(pg))]
int <- int[,-grep("Blank", names(int))]

int[int == 0] <- NA

# remove all 0 columns
int <- int[,colSums(is.na(int)) < nrow(int)]

# number of observations
num.obs <- apply(int, 1, function(x) sum(!is.na(x)))

median <- apply(int, 1, function(x) median(x, na.rm = T))

# split proteins into deciles by median intensity
decile <- ntile(-median, 10)

# bind df
df <- cbind.data.frame(pg$Sequence, num.obs, median, decile)

# allCounts
count <- ggplot(df, aes(num.obs))+
  geom_histogram(binwidth = 1, color = "black", fill = "white")+
  scale_x_reverse()+
  labs(x = NULL, y = "Peptide Count")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20))
  
#svg(file = paste(dir2, "ProteinAbundanceCurve_bottomPlot.svg", sep = ""), height = 5, width = 7)
scurve <- ggplot(df, aes(num.obs, median, color = factor(decile)))+
  geom_point(shape = 1, alpha = 0.6, size = 1)+
  scale_y_log10()+
  scale_color_brewer(palette = "Spectral")+
  scale_x_reverse()+
  labs(x = "Number of observations", y = "Median Intensity", color = "Decile")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))
#dev.off()


svg(file = "R_figures/ProteinAbundanceCurve_legend_peptide.svg")
ggarrange(scurve, count, nrow = 2, align = "hv", heights = c(1, 0.5), common.legend = TRUE)
dev.off()

rm(list = ls())
q()
