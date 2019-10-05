#!/usr/local/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

print(args)

for (arg in args) {
  data <- read.table(arg, sep = "\t", header = T)
  
  print(head(data[,1:20]))
}