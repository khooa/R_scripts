source("D:/R_scripts/mqlib_prot.R")
theme_set(plot_theme())

exo <- pgRead("D:/PCa_Discovery/Data_analysis_files/txt_190917_mbr_exo/R_tableOutput/proteinCounts_mbr.txt")
exo$sample <- gsub("B1297", "B1287", exo$sample)

mv <- pgRead("D:/PCa_Discovery/Data_analysis_files/txt_190925_mv_mbr//R_tableOutput/proteinCounts_mbr.txt")

exo_id <- getSampleIDs(exo$sample)
exo$sample <- str_glue_data(exo_id, "{DRE}_{patient}_{timepoint}")

mv_id <- getSampleIDs(mv$sample)
mv$sample <- str_glue_data(mv_id, "{DRE}_{patient}_{timepoint}")

names(exo) <- gsub("value", "exo", names(exo))
names(mv) <- gsub("value", "mv", names(mv))

exo <- subset(exo, variable == "MS2")
mv <- subset(mv, variable == "MS2")

new_df <- merge(exo, mv, by = "sample", all = TRUE)

new_df$variable.x <- NULL
new_df$variable.y <- NULL

rsq <- summary(lm(new_df$exo ~ new_df$mv))$r.squared

ggplot(new_df, aes(exo, mv))+
  geom_point()+
  geom_smooth(method = "lm", formula = y~x, se = FALSE)
