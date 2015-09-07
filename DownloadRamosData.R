#download the Ramos et al GEO dataset

setwd("~/Testing123/PairedOAAnalysis")


library(GEOquery)

gds <- getGEO("GSE57218")

gset <-gds$GSE57218_series_matrix.txt.gz

save(gset,file="GSE57218.RData")