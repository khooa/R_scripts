#!/usr/bin/env Rscript

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

write.table(mdf, "R_tableOutput/proteinCounts_mbr.txt", sep = "\t", row.names = F)

# sample order
sample_order <- subset(mdf, variable == "MS2")
sample_order <- sample_order[order(sample_order$value),]

svg("R_figures/protCount_barplot_lfq_withMatching.svg", width = 20, height = 4)
ggplot(mdf, aes(sample, value, fill = variable))+
  geom_bar(stat = "identity", position = position_dodge())+
  scale_x_discrete(limits = sample_order$sample, labels = sample_order$sample)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))+
  geom_text(aes(label = value), position = position_dodge(width = 0.9), angle = 90, hjust = 1.2, vjust = 0.3, size = 2)+
  scale_fill_manual(values = c("#4682B4", "#a9a9a9"))+
  labs(x = NULL, y = "Number of Proteins", fill = "ID Type")
dev.off()

print("protCount_barplot_lfq_withMatching.svg done")

svg("R_figures/protCount_barplot_stacked_lfq_withMatching.svg", width = 20, height = 4)
ggplot(mdf, aes(sample, value, fill = variable))+
  geom_bar(stat = "identity", position = position_stack(rev = T))+
  scale_x_discrete(limits = sample_order$sample, labels = sample_order$sample)+
  geom_text(aes(label = value), position = position_stack(rev = T), angle = 90, hjust = 1.2, vjust = 0.3, size = 2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))+
  scale_fill_manual(values = c("#4682B4", "#a9a9a9"))+
  labs(x = NULL, y = "Number of Proteins", fill = "ID Type")
dev.off()

print("protCount_barplot_stacked_lfq_withMatching.svg done")
