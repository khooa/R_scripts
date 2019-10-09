#!/usr/bin/env Rscript

source("D:/R_scripts/mqlib_prot.R")
theme_set(plot_theme())

# create folders
dir.create("R_tableOutput")
dir.create("R_figures")

##-----------MS PERFORMANCE -----------------##

summary <- read.table("summary.txt", sep = "\t", header = T)

summary <- filter(summary, Experiment != "")

summary <- summary[,c("Experiment", "MS", "MS.MS", "MS.MS.Identified", "MS.MS.Identified....", 
                      "Peaks", "Peaks.Sequenced", "Peaks.Repeatedly.Sequenced")] # select columns

# split sample name into fraction, sample, run, lot and gleason.score
sample.id <- getSampleIDs(summary$Experiment)

# bind dataframes
summary <- cbind.data.frame(sample.id, summary)

##------------PROTEIN GROUPS----------------##
pg <- read.table("proteinGroups.txt", sep = "\t", header = T)

#filter data
pg_f <- pgRemoveCon(pg)
pg_f <- pgRemoveRev(pg_f)
pg_f <- pgRemoveIDbySite(pg_f)
pg_f <- pgMinPept(pg_f)


# get gene names, uniprot ids
protein.id <- getInfo(pg_f)

# Merge info df with intensity
makeDF <- function(x) {
  x[x == 0] <- NA

  median <- apply(x, 1, function(y) median(y, na.rm = T))

  num.obs <- apply(x, 1, function(x) sum(!is.na(x)))

  df <- cbind.data.frame(protein.id, median, num.obs, x)

  return(df)
}

# get ibaq data
ibaq <- pgGetiBAQ(pg_f)

names(ibaq) <- paste("ibaq", names(ibaq), sep = "_")

ibaq_df <- makeDF(ibaq)

# get lfq data
lfq <- pgGetLFQ(pg_f)

names(lfq) <- paste("lfq", names(lfq), sep = "_")

lfq_df <- makeDF(lfq)

##------------PEPTIDES----------------##
# read table
pep <- pgRead("peptides.txt")

# filter data
pep <- pgRemoveRev(pep)
pep <- pgRemoveCon(pep)
pep.lfq <- pgGetLFQ(pep)

names(pep.lfq) <- paste("lfq", names(pep.lfq), sep = "_")

peptide.info <- pep[,c("Sequence", "Missed.cleavages", "Charges", "Score", "PEP", "Leading.razor.protein", "Gene.names")]

pep_lfq_df <- cbind.data.frame(peptide.info, pep.lfq)

##-------PROTEIN and PEPTIDE COUNTS ------------##
##------- based on imputed values---------------##

pg.count <- apply(lfq, 2, function(x) length(x[x>0]))
pep.count <- apply(pep.lfq, 2, function(x) length(x[x>0]))

count <- cbind.data.frame(names(lfq), pg.count, pep.count)
names(count) <- c("Experiment", "prot.count", "pep.count")

summary <- merge(summary, count, by = "Experiment")

##-----WRITE DATA TABLES------------##
tableList <- list(summary, ibaq_df, lfq_df, pep_lfq_df)
names(tableList) <- c("summary", "ibaq_df", "lfq_df", "pep_lfq_df")

sapply(names(tableList),
       function(x) write.table(tableList[[x]], file = paste("R_tableOutput/", x, ".txt", sep = ""), sep = "\t", col.names = T, row.names = F))
