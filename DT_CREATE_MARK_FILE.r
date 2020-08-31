library(purrr)
library(dplyr)
library(org.Hs.eg.db)
library(annotate)

setwd('f:/mmrf_data')

df <- read.csv('IGV_downloads/MMRF_CoMMpass_IA11a_IGV_CNA_Exome.seg', sep="\t", header=TRUE, skip = 1)

head(df)

rm(dfFinal)

dfFinalS <- df %>% dplyr::select(ID, chrom, loc.start)

dfFinalE <- df %>% dplyr::select(ID, chrom, loc.end)

names(dfFinalS) <- c("Marker Name", "Chrom", "Pos")
names(dfFinalE) <- c("Marker Name", "Chrom", "Pos")

head(dfFinalS)
head(dfFinalE)

dfFinal <- bind_rows(dfFinalS, dfFinalE)

nrow(dfFinalS)

nrow(dfFinalE)

nrow(dfFinal)

head(dfFinal,20)

dfFinal <- dfFinal %>% arrange(`Marker Name`, Chrom, Pos)

head(dfFinal,40)


outputFile <- "DT_R/cbio_files/mark_file.txt"

file.remove(outputFile)
file.create(outputFile)

write.table(dfFinal, file=outputFile, sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, append = TRUE, na="") 

