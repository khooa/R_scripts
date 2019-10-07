library("ggplot2")
library("reshape2")
library("extrafont")
library("circlize")
library("data.table")
library("scales")
library("dplyr")
library("gdata")
library("limma")
library("ggrepel")
library("RColorBrewer")
library("ggpubr")
library("VennDiagram")
library("ComplexHeatmap")
library("UpSetR")
library("stringr")

#'------------------------------------------------------------------------------------------
#' Colours
#' -----------------------------------------------------------------------------------------

# colours
cell_colours <- c("LNCaP" = "#e41a1c",
                  "22Rv1" = "#ff7f00",
                  "DU145" = "#1f78b4",
                  "PC3" = "#33a02c",
                  "RWPE1" = "#2f4f4f",
                  "NA" = "#a9a9a9")

ev_colours <- c("cell" = "#b3b3b3",
                "ev" = "#e41a1c",
                "non_ev" = "#377eb8",
                "NA" = "#a9a9a9")

tissue_colours <- c("prostate" = "#b22222",
                    "kidney" = "#f0e68c",
                    "urinary_tract" = "#a6761d",
                    "ovary" = "#66a61e",
                    "tongue" = "#7570b3",
                    "NA" = "#a9a9a9")
#'------Function to extract colours

tissue_cols <- function(...) {
  cols <- c(...)
  
  if(is.null(cols))
    return(tissue_colours)
  
  tissue_colours[cols]
}

cell_cols <- function(...) {
  cols <- c(...)
  
  if(is.null(cols))
    return(cell_colours)
  
  cell_colours[cols]
}

ev_cols <- function(...) {
  cols <- c(...)
  
  if(is.null(cols))
    return(ev_colours)
  
  ev_colours[cols]
}

#'-----------------------------------------------------------------------------------------
#' Marker proteins------------------------------------------------------------------------
#' ----------------------------------------------------------------------------------------

prostate_markers <- c("KLK3", "FOLH1", "KLK2", "MSMB", "ACPP", "NKX3-1", "NEFH", "AR", "RLN1", "KLK4", "SLC45A3", "STEAP2", "TGM4", "SP8")
ev_markers <- c("CD9", "CD81", "CD63", "TSG101", "ARRDC1", "ANXA2", "ANXA5", "ANXA11", "PDCD6IP", "SDCBP", "FLOT1", "FLOT2",
                "ITGB2", "FN1")
nonev_markers <- c("UMOD", "APOA1", "LMNA", "CANX", "GAPDH", "ALB", "CYC1", "GOLGA2", "HSP90B1", "CALR")

#'------------------------------------------------------------------------------------------
#' Plot theme
#' -----------------------------------------------------------------------------------------

plot_theme <- function() {
  theme_classic() %+replace%
    theme(
    axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = "bold"),
        title = element_text(size = 18, face = "bold"),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18, face = "bold")
        )
}

#'-------------------------------------------------------------------------------------------
#' Read ProteinGroup.txt file
#'-------------------------------------------------------------------------------------------
pgRead <- function(file) {
  pg <- data.frame(fread(file, sep="\t", header=T, stringsAsFactors=FALSE, integer64 = "numeric", data.table = FALSE))
  rownames(pg) <- pg$id
  return(pg)
} 

#'-------------------------------------------------------------------------------------------
#' Remove reversed protein
#'-------------------------------------------------------------------------------------------
pgRemoveRev <- function(pg) {
  if (sum(!is.na(pg$Reverse))==0 | sum(pg$Reverse=="+")==0) {
    print("No revesed protein to be removed ... ")
  } else {
    pg <- pg[(!pg$Reverse=="+"),]
  }
  return(pg)
}

#'-------------------------------------------------------------------------------------------
#' Remove potential contaminant
#'-------------------------------------------------------------------------------------------
pgRemoveCon <- function(pg) {
  if (sum(!is.na(pg$Potential.contaminant))==0 | sum(pg$Potential.contaminant=="+")==0) {
    print("No contaminant protein to be removed ... ")
  } else {
    pg <- pg[(!pg$Potential.contaminant=="+"),]
  }
  return(pg)
}

#'-------------------------------------------------------------------------------------------
#' Remove only identified by site
#'-------------------------------------------------------------------------------------------
pgRemoveIDbySite <- function(pg) {
  if (sum(!is.na(pg$Only.identified.by.site))==0 | sum(pg$Only.identified.by.site=="+")==0) {
    print("No Only.identified.by.site protein to be removed ... ")
  } else {
    pg <- pg[(!pg$Only.identified.by.site=="+"),]
  }
  return(pg)
}

#'-------------------------------------------------------------------------------------------
#' Filter by minimum number of peptides
#'-------------------------------------------------------------------------------------------
pgMinPept <- function(pg, n=2){
  ru <- pg[grep("Razor...unique.peptides.", colnames(pg))]
  pg <- pg[apply(ru>=n, 1, sum)>0,]
  return(pg)
}

#'-------------------------------------------------------------------------------------------
#' Get columns with LFQ Intensity from ProteinGroup dataframe
#'-------------------------------------------------------------------------------------------
pgGetLFQ <- function(pg) {
  lfq <- pg[grep("LFQ.intensity.", colnames(pg))]
  colnames(lfq) <- gsub("LFQ.intensity.", "", colnames(lfq)) 
  return(lfq)
}

#'-------------------------------------------------------------------------------------------
#' Get columns with iBAQ from ProteinGroup dataframe
#'-------------------------------------------------------------------------------------------
pgGetiBAQ <- function(pg) {
  ibaq <- pg[grep("iBAQ.", colnames(pg))]
  colnames(ibaq) <- gsub("iBAQ.", "", colnames(ibaq)) 
  return(ibaq)
}

#'-------------------------------------------------------------------------------------------
#' Impute median adjusted iBAQ to LFQ
#'-------------------------------------------------------------------------------------------
imputeiBAQtoLFQ <- function(data, logTrans=F, normMethod="median", rowStatus=F) {
  lfq.col <- grep("^LFQ.intensity.", colnames(data))
  ibaq.col <- grep("^iBAQ.", colnames(data))
  lfq <- data[,lfq.col]
  ibaq <- data[,ibaq.col]
  lfq[lfq==0]<-NA
  ibaq[ibaq==0]<-NA
  lfq <- log2(lfq)
  ibaq <- log2(ibaq)
  lfq[is.na(lfq)] <- 0
  ibaq[is.na(ibaq)] <- 0
  lfq.s <- apply((lfq>0)*1, 1, sum)
  ibaq.s <- apply((ibaq>0)*1, 1, sum)
  lfq.complite <- lfq[lfq.s==length(lfq.col),]
  ibaq.complite <- ibaq[lfq.s==length(lfq.col),]
  if (normMethod=="median") {
    lfq.median <- apply(lfq.complite, 2, median)
    ibaq.median <- apply(ibaq.complite, 2, median)
    m.dif <- lfq.median - ibaq.median 
  } else if(normMethod=="mean"){
    lfq.mean <- apply(lfq.complite, 2, mean)
    ibaq.mean <- apply(ibaq.complite, 2, mean)
    m.dif <- lfq.mean - ibaq.mean 
  } else {
    m.dif <- 0
  }
  ibaq[ibaq==0]<-NA
  for(i in 1:length(ibaq.col)) {
    ibaq[,i] <- ibaq[,i] + m.dif[i]
  }
  ibaq[is.na(ibaq)] <- 0
  lfq.adj.ibaq<-lfq
  rowsToImpute<-apply(((lfq>0)*1),1,sum)<dim(lfq)[2]
  lfq.adj.ibaq[rowsToImpute,]<-ibaq[rowsToImpute,]
  if (logTrans==F) {
    lfq.adj.ibaq <- data.frame(2^lfq.adj.ibaq) 
    lfq.adj.ibaq[lfq.adj.ibaq==1] <- 0
  }
  if (rowStatus==T) {
    lfq.adj.ibaq$rowStatus<-"lfq"
    lfq.adj.ibaq$rowStatus[rowsToImpute]<-"ibaq"
  }
  #pIDs <- grep("^Protein.IDs", colnames(data))
  #rownames(lfq.adj.ibaq)<-data[,pIDs]
  rownames(lfq.adj.ibaq)<-data$id
  names(lfq.adj.ibaq) <- gsub("LFQ.intensity.", "LFQ.adj.iBAQ.", names(lfq.adj.ibaq))
  return(lfq.adj.ibaq)
}

#'-------------------------------------------------------------------------------------------
#' Impute from lower distribution (Perseus method)
#'-------------------------------------------------------------------------------------------
imputeLD <- function(data, shift=1.8, size=0.2) {
  set.seed(8)
  data[data==0]<-NA
  mdf <- melt(data, na.rm=T)
  n<-round(dim(mdf)[1])
  sdev <- sd(mdf$value)
  m <- mean(mdf$value) - (sdev * shift)
  d<-rnorm(n*size, mean=m, sd=sdev*size)
  data[is.na(data)] <- d[round(runif(sum(is.na(data)),1,n*size))]
 return(data)
}

#'-------------------------------------------------------------------------------------------
#' Transform data either Log2 or Log10
#'-------------------------------------------------------------------------------------------
transform <- function(data, base="log2", na=F) {
  data[data==0] <- NA
  if (base=="log2") {
    data <- log2(data)
  } else if (base=="log10") {
    data <- log10(data)
  }
  if (na==F) {
    data[is.na(data)]<-0
  }
  return(data)
}

#'------------------------------------
#'-- bind vectors of different lengths
#'------------------------------------
cfun <- function(L) {
  pad.na <- function(x,len) {
    c(x,rep(NA,len-length(x)))
  }
  maxlen <- max(sapply(L,length))
  do.call(data.frame,lapply(L,pad.na,len=maxlen))
}

#'------------------------------------
#'-- Get ID--------------------------
#'------------------------------------
getSampleIDs <- function(x) {
  
  id <- as.data.frame(str_split_fixed(x, "_", 4))
  names(id) <- c("DRE", "EV", "patient", "timepoint")
  
  id$DRE <- gsub("DRE", "", id$DRE)
  
  return(id)
}

getInfo <- function(x) {
  genes <- as.data.frame(str_split_fixed(x$Gene.names, ";", 2))$V1
  main.isoform <- as.data.frame(str_split_fixed(x$Majority.protein.IDs, ";", 2))$V1
  
  minor.isoforms <- as.data.frame(str_split_fixed(x$Majority.protein.IDs, ";", 2))$V2
  
  protein <- as.data.frame(str_split_fixed(main.isoform, "-", 2))$V1
  
  protein.names <- x$Protein.names
  
  rzr.uniq.pep <- grep("^Razor...unique.peptides", names(x))
  
  unique.peptides <- x[,rzr.uniq.pep]
  
  id <- cbind.data.frame(genes, protein, main.isoform, minor.isoforms, protein.names, unique.peptides)
  
  return(id)
}
