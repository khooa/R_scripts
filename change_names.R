library("data.table")

setwd("D:/PCa_Discovery/Data_analysis_files/txt_exo_16Jan2019/")

pg <- data.frame(fread("proteinGroups.txt", sep="\t", header=T, stringsAsFactors=FALSE, integer64 = "numeric", data.table = FALSE))
ev <- data.frame(fread("D:/PCa_Discovery/Data_analysis_files/Patient_info/EV_clinicalCovariates.txt", sep="\t", header=T, stringsAsFactors=FALSE, integer64 = "numeric", data.table = FALSE))

pg.names <- colnames(pg[,grep("Intensity.", colnames(pg))])
pg.names <- data.frame(old=sub("Intensity.", "", pg.names), new=NA, E=NA, R=NA, P=NA, G1=NA, G2=NA, AS=NA)

for (i in 1:length(pg.names$old)) {
  parts <- strsplit(as.vector(pg.names$old[i]), "_")[[1]]
  pg.names$E[i] <- parts[1]
  pg.names$R[i] <- parts[2]
  pg.names$P[i] <- parts[3]
  pg.names$G1[i] <- parts[4]
  pg.names$G2[i] <- parts[5]
  pg.names$AS[i] <- parts[6]
}

for (i in 1:length(pg.names$old)) {
  if (nchar(pg.names$P[i])<5) {
    pg.names$P[i] <- gsub("P","P0",pg.names$P[i])
  }
}

pg.names <- merge(pg.names, ev[,c("ISUP","fraction")], by.x = "E", by.y = "fraction")
for (i in 1:length(pg.names$old)) {
  if (is.na(pg.names$AS[i])) {
    pg.names$new[i] <- paste(c(pg.names$E[i], "_", pg.names$P[i], "_ISUP", pg.names$ISUP[i]), collapse = "")
  } else {
    pg.names$new[i] <- paste(c(pg.names$E[i], "_", pg.names$P[i], "_ISUP",pg.names$ISUP[i], "_", pg.names$AS[i]), collapse = "")
  }
}

# replace old sample names to new in protein group file
for (i in 1:length(pg.names$old)) {
  colnames(pg) <- gsub(pg.names$old[i], pg.names$new[i], colnames(pg))
}

fwrite(pg, "R_tableOutput/proteinGroups.txt", sep = "\t")
