source("D:/R_scripts/mqlib_prot.R")

#setwd("D:/PCa_Discovery/Data_analysis_files/subset_prostate_proteins/")

#'--------------------------------
#'---------Functions--------------
#'--------------------------------

getGenes <- function(x) {
  
  data <- pgRead(x)
  genes <- data$genes
  
  return(genes)
  
}

prostateGeneDF <- function(x, dataset) {
  gene_df <- cbind.data.frame(x, rep(dataset, length(x)))
  names(gene_df) <- c("genes", "dataset")
  
  return(gene_df)
}

#'-------------------------------------------------
#'--------Read tk lab files & get proteins --------
#'-------------------------------------------------

# read in files
filelist <- list.files(pattern = "pca*")

genelist_data <- lapply(filelist, function(x) getGenes(x))

# give names to genes in list

names(genelist_data) <- gsub(".txt", "", filelist)

# add descriptors

genelist_prot_df <- list()

for (i in seq_along(genelist_data)) {
  
  create_df <- data.frame(genes = genelist_data[i], dataset = rep(names(genelist_data)[i], length(genelist_data[i])))
  
  names(create_df) <- c("genes", "dataset")
  
  genelist_prot_df[[length(genelist_prot_df) + 1]] <- create_df
  
}

names(genelist_prot_df) <- names(genelist_data)

#'-------------------------------------------------
#'--------Get HPA prostate proteins ---------------
#'-------------------------------------------------

hpa <- pgRead("hpa_prostate.tsv")
hpa[hpa == ""] <- NA

attach(hpa)

hpa_cancer_rna <- subset(hpa, RNA.cancer.category != "Not detected" | !is.na(RNA.cancer.category))
hpa_cellLine_rna <- subset(hpa, RNA.cell.line.category != "Not detected" | !is.na(RNA.cell.line.category))
hpa_tissue_rna <- subset(hpa, RNA.tissue.category != "Not detected" | !is.na(RNA.tissue.category))
hpa_tissue_protein <- subset(hpa, Evidence == "Evidence at protein level")

detach(hpa)

hpalist <- list(hpa_cancer_rna, hpa_cellLine_rna, hpa_tissue_protein, hpa_tissue_rna)

genelist_hpa <- lapply(hpalist, function(x) x$Gene)
names(genelist_hpa) <- c("hpa_cancer_rna", "hpa_cellLine_rna", "hpa_tissue_protein", "hpa_tissue_rna")

genelist_hpa_df <- list()

for (i in seq_along(genelist_hpa)) {
  
  # i = 2
  
  create_df <- data.frame(genes = genelist_hpa[i], dataset = rep(names(genelist_hpa)[i], length(genelist_hpa[i])))
  
  names(create_df) <- c("genes", "dataset")
  
  genelist_hpa_df[[length(genelist_hpa_df) + 1]] <- create_df
  
}

names(genelist_hpa_df) <- names(genelist_hpa)

#'--------------------------------------------------
#'----- Get CCLE prostate genes -----------------
#'--------------------------------------------------

ccle <- pgRead("ccle_rsem_prostate_genes.txt")

genes_ccle <- prostateGeneDF(ccle$genes, "ccle")

#'--------------------------------------------------
#'----- Get GTEX prostate genes -----------------
#'--------------------------------------------------

gtex <- pgRead("GTEx_v8_gene_median_tpm.gct")

gtex[gtex == 0] <- NA

gtex_prostate <- subset(gtex, !is.na(Prostate))

genes_gtex <- prostateGeneDF(gtex$Description, "gtex")

#'--------------------------------------------------
#'----- Get Gato et al prostate proteins -----------
#'--------------------------------------------------

gato <- pgRead("pca_gato_tissue.txt")



#'------------------------------------------------------------------
#'----- Make dataframe of genes detected in each dataset ------------
#'------------------------------------------------------------------

# Merge lists of genes

dataset_list <-c(genelist_prot_df, genelist_hpa_df, list(genes_ccle, genes_gtex))

names(dataset_list) <- c(names(genelist_prot_df), names(genelist_hpa_df), "ccle", "gtex")

# Merge dataframes by gene name

dataset_df <- Reduce(function(x, y) merge(x, y, all = TRUE), dataset_list)

# remove rows with no gene names
dataset_df <- subset(dataset_df, genes != "")

# reformat into wide format
dataset_df_wide <- dcast(dataset_df, genes ~ dataset)

dataset_df_wide[-1][dataset_df_wide[-1] > 0] <- 1

# write out file

write.table(dataset_df_wide, file = "D:/PCa_Discovery/Data_analysis_files/subset_prostate_proteins/R_tableOutput/prostate_proteins_unfiltered.txt",
            sep = "\t", row.names = F)
