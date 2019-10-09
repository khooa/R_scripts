setwd("D:/PCa_Discovery/Data_analysis_files/CellLine/txt_all_cellLine_EV/")

library(data.table)
library(dplyr)
library(stringr)
library(reshape2)
library(VennDiagram)
library(UpSetR)

plot_theme <- theme(axis.text = element_text(size = 14),
                    axis.title = element_text(size = 14, face = "bold"),
                    title = element_text(size = 18, face = "bold"),
                    strip.text.x = element_text(size = 14),
                    strip.text.y = element_text(size = 14),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 18, face = "bold"))

# colours
cell_colours <- c("LNCaP" = "#e41a1c",
                  "22Rv1" = "#ff7f00",
                  "DU145" = "#1f78b4",
                  "PC3" = "#33a02c",
                  "RWPE1" = "#000000")

# read in files

ccle <- read.table("R_tableOutput/10_CCLE_all_prostate_proteins.txt", sep = "\t", header = T)
ccle$CCLE.category <- ifelse(ccle$num_enriched == 25, "prostate_enriched",
                             ifelse(ccle$num_enriched < 25 & ccle$num_enriched >= 20, "prostate_enhanced", "expressed_in_all"))
ccle <- ccle[,c("genes", "CCLE.category")]

hpa <- read.table("R_tableOutput/11_HPA_prostate.txt", sep = "\t", header = T)

cpcg_tissue <- read.table("D:/PCa_Discovery/Data_analysis_files/datasets/CPCGENE_PCa_Tissue/R_tableOutput/pg_imp_df.txt", sep = "\t", header = T)
cpcg_tissue <- as.character(cpcg_tissue$genes)

pca_wcl <- read.table("D:/PCa_Discovery/Data_analysis_files/CellLine/txt_190918_pca_wcl/R_tableOutput/pg_imp_df.txt", sep = "\t", header = T)
pca_wcl <- as.character(pca_wcl$genes)

listInput <- list(CCLE = ccle$genes, HPA = hpa$Gene, CPCG_PCa_tissue = cpcg_tissue, pca_wcl = pca_wcl)

svg("R_figures/13_UpsetPlot_shared_prostate_proteins.svg", height = 6, width = 8)
upset(fromList(listInput), order.by = "freq")
dev.off()

overlap <- calculate.overlap(listInput)

shared_prot <- overlap$a6

#write.table(shared_prot, "R_tableOutput/13_shared_prostate_proteins.txt", sep = "\t", row.names = F)
