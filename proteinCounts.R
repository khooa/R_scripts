#setwd("D:/PCa_Discovery/Data_analysis_files/txt_190917_mbr_exo/")

source("D:/PCa_Discovery/R_scripts/mqlib_prot.R")
theme_set(plot_theme())

pg <- read.table("proteinGroups.txt", sep = "\t", header = T)

pg_f <- pgRemoveCon(pg)
pg_f <- pgRemoveRev(pg_f)
pg_f <- pgRemoveIDbySite(pg_f)
pg_f <- pgMinPept(pg_f)

mbr <- grep("^Identification.type.", names(pg_f))

pg_id <- pg_f[,mbr]

by_MS2 <- apply(pg_id, 2, function(x) sum(x == "By MS/MS"))

by_matching <- apply(pg_id, 2, function(x) sum(x == "By matching"))

df <- cbind.data.frame(names(pg_id), by_MS2, by_matching)
names(df) <- c("sample", "MS2", "matching")

mdf <- melt(df, id.vars = "sample")

mdf$sample <- gsub("Identification.type.", "", mdf$sample)

# split id vars
id_var <- getSampleIDs(mdf$sample)
id_var$group <- ifelse(str_detect(id_var$patient, "^B"), "benign", "cancer")

mdf_id <- cbind.data.frame(id_var, mdf)

# sample order
sample_order <- subset(mdf_id, variable == "MS2")
sample_order <- sample_order[order(sample_order$value),]

svg("R_figures/protCount_barplot_lfq_withMatching.svg", width = 16, height = 4)
ggplot(mdf_id, aes(sample, value, fill = variable))+
  geom_bar(stat = "identity", position = position_dodge())+
  scale_x_discrete(limits = sample_order$sample, labels = mdf_id$patient)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))
dev.off()

svg("R_figures/protCount_barplot_stacked_lfq_withMatching.svg", width = 16, height = 4)
ggplot(mdf_id, aes(sample, value, fill = variable))+
  geom_bar(stat = "identity", position = position_stack(rev = T))+
  scale_x_discrete(limits = sample_order$sample, labels = mdf_id$patient)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))
dev.off()