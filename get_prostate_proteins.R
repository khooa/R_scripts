source("D:/R_scripts/mqlib_prot.R")
theme_set(plot_theme())

#setwd("D:/PCa_Discovery/Data_analysis_files/subset_prostate_proteins/")

data <- pgRead("R_tableOutput/prostate_proteins_unfiltered.txt")

#'----------------------------------------------------------------------------------
#'-------Get set of common prostate proteins from whole cell/tissue data------------
#'-----------------------------------------------------------------------------------

data_subset <- data[,c("genes", "pca_cellLine_wcl", "pca_cpcg_tissue", "hpa_cancer_rna", "hpa_cellLine_rna",
                       "hpa_tissue_protein", "hpa_tissue_rna", "ccle", "gtex")]

# Count number of datasets each gene is detected in

data_subset$sum <- apply(data_subset[-1], 1, function(x) sum(x))

pdf("R_figures/histogram_dataset_count.pdf", height = 5, width = 5)
hist(data_subset$sum,
     main = "Number of datasets each gene is detected in",
     xlab = "Number of datasets",
     ylab = "Number of genes",
     breaks = seq(0,8,1))
dev.off()

# Get genes in more than 5 datasets

data_prostate_filtered <- subset(data_subset, sum > 4)

# write table to file
write.table(data_prostate_filtered, file = "R_tableOutput/prostate_proteins_filtered_191009.txt", sep = "\t", row.names = F)

#---------------------------------------------------------------------------
#----- Make plot of "likely prostate proteins" vs "prostate proteins"-------
#---------------------------------------------------------------------------

prostate_prot_count <- data.frame(matrix(ncol = 3, nrow = 8))
names(prostate_prot_count) <- c("dataset", "all datasets (8)", ">4 datasets")

prostate_prot_count[1] <- names(data_subset[,-c(1, 10)])
prostate_prot_count[2] <- apply(data_subset[,-c(1, 10)], 2, function(x) sum(x))
prostate_prot_count[3] <- apply(data_prostate_filtered[,-c(1, 10)], 2, function(x) sum(x))

prostate_prot_count[2] <- prostate_prot_count[3] - prostate_prot_count[2]

m_prostate_prot_count <- melt(prostate_prot_count, id.vars = "dataset")

svg("R_figures/barplot_genes_filtered.svg", height = 5, width = 8)
ggplot(m_prostate_prot_count, aes(dataset, value, fill = variable))+
  geom_bar(stat = "identity", position = position_stack())+
  scale_fill_manual(values = c("#a9a9a9", "#377eb8"))+
  labs(x = "Number of genes", y = NULL, fill = NULL)+
  geom_text(aes(label = value), vjust = 1, hjust = 0.5, fontface = "bold")+
  coord_flip()+
  theme(legend.position = "top")
dev.off()

#---------------------------------------------------------------------------
#----- Make plot of "likely prostate proteins" in proteomics data ----------
#--------------------------------------------------------------------------

proteomics_data <- data[,c(1, grep("^pca", names(data)))]
proteomics_data_f <- proteomics_data[proteomics_data$genes %in% data_prostate_filtered$genes,]

proteomics_prot_count <- data.frame(matrix(ncol = 3, nrow = 7))
names(proteomics_prot_count) <- c("dataset", "total_prot", "prostate_prot")

proteomics_prot_count[1] <- names(proteomics_data_f[-1])
proteomics_prot_count[2] <- apply(proteomics_data[-1], 2, function(x) sum(x))
proteomics_prot_count[3] <- apply(proteomics_data_f[-1], 2, function(x) sum(x))

proteomics_prot_count$percent_difference <- ((proteomics_prot_count$total_prot - proteomics_prot_count$prostate_prot)/proteomics_prot_count$total_prot)*100

m_proteomics_count <- melt(proteomics_prot_count, id.vars = "dataset")

# make filled barplot

xaxis_order <- proteomics_prot_count[order(-proteomics_prot_count$percent_difference),]

ggplot(subset(m_proteomics_count, variable != "total_prot"), aes(dataset, value, fill = variable))+
  geom_bar(stat = "identity", position = position_fill(reverse = T))+
  labs(x = NULL, y = "Percentage of proteins", fill = NULL)+
  scale_fill_manual(values = c("#377eb8", "#a9a9a9"))+
  scale_x_discrete(limits = xaxis_order$dataset)+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
