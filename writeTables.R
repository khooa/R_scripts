source("D:/PCa_Discovery/R_scripts/mqlib_prot.R")
theme_set(plot_theme())

# create folders
dir.create("R_tableOutput")
dir.create("R_figures")

##-----------MS PERFORMANCE -----------------##

# Read in summary.txt file

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

# filter data
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
ibaq_df <- makeDF(ibaq)

# get lfq data
lfq <- pgGetLFQ(pg_f)
lfq_df <- makeDF(lfq)

##------------PEPTIDES----------------##
# read table
pep <- pgRead("peptides.txt")

# filter data
pep <- pgRemoveRev(pep)
pep <- pgRemoveCon(pep)
pep.lfq <- pgGetLFQ(pep)

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
#'------------------------------------------------------------
#'--------------- join with clinical covariates---------------
#'------------------------------------------------------------

info <- read.table("D:/PCa_Discovery/Data_analysis_files/Patient_info/EV_patient_info.txt", sep = "\t", header = T)

# group ages
info$age.group <- cut(info$Age.Dx, breaks = c(-Inf, 49, 59, 69, Inf), labels = c("<50", "50-59", "60-69", ">= 70"))

info$Tcat <- ifelse(str_detect(info$cT, "T1"), "T1",
                    ifelse(str_detect(info$cT, "T2"), "T2",
                           ifelse(str_detect(info$cT, "T3"), "T3",
                                  ifelse(str_detect(info$cT, "T4"), "T4", NA))))

info$psa <- cut(info$PSA.Dx, breaks = c(-Inf, 9, 20, Inf), labels = c("<10", "10-20", ">20"))

info <- filter(info, EV.subtype == "EXO")

info$cT <- gsub("T3b;T4", "T4", info$cT)

# load pg file
pg <- read.table("R_tableOutput/lfq_df.txt", sep = "\t", header = T)

# sort by median intensity and remove duplicate genes
pg <- pg[order(pg$median),]

pg <- pg[!duplicated(pg$genes),]

# get intensities

int.cols <- c(grep("^postDRE_", names(pg)), grep("^preDRE_", names(pg)))

pg.int <- pg[,int.cols]

rownames(pg.int) <- pg$genes

# transpose data for PCA plot (colnames = genes, rownames = samples)
t.pg.int <- as.data.frame(t(pg.int))

# add Experiment name as identifier
t.pg.int$experiment <- rownames(t.pg.int)

#'--------------------------------------------------------------------------------
#'------ Merge samples with batch data--------------------------------------------
#'--------------------------------------------------------------------------------

samples <- data.frame(experiment = rownames(t.pg.int), DRE=NA, EV=NA, patient=NA, timepoint=NA) # create empty dataframe with colnames

for (i in 1:length(rownames(t.pg.int))) {
  parts <- strsplit(as.vector(rownames(t.pg.int)[i]), "_")[[1]]
  samples$DRE[i] <- gsub("DRE", "", parts[1])
  samples$EV[i] <- parts[2]
  samples$patient[i] <- parts[3]
  samples$timepoint[i] <- parts[4]
}

# add "experiment" column to info df
info$experiment <- str_glue_data(info, "{DRE}DRE_EXO_{PID}_{Longitudinal}")

# Merge covariates and intensity table
data <- merge(info, t.pg.int, by = "experiment")

# write table with covariates
write.table(data, "R_tableOutput/pg_imp_covariates.txt", sep = "\t", row.names = F)

rm(list=ls())