#!/usr/bin/env Rscript
library(rtracklayer)
gff <- import(commandArgs(T)[1], format='GFF')
tabular <- as.data.frame(gff)
head(tabular)
# some columns are lists, we need to convert them to vectors  before writing to file
for (i in 1:ncol(tabular)){
  if (is.list(tabular[[i]])){
    tabular[[i]] <- sapply(tabular[[i]], function(x) paste(x, collapse = ";"))
  }
}
write.table(tabular, file = commandArgs(T)[2], quote=FALSE, sep="\t", row.names=FALSE)

