library(purrr)
library(dplyr)
library(org.Hs.eg.db)
library(annotate)

setwd('C:/Users/abhmalat/OneDrive - Indiana University/cbio_MM/')

# Fix quotes for height
df <- read.csv('data_timeline_lab_test_CMG_MM.txt', sep="\t", header=TRUE, quote = "")
write.table(df, file='data_timeline_lab_test_CMG_MM_formatted.txt', sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")


df_treatment <- read.csv('data_timeline_treatment_CMG_MM.txt', sep="\t", header=TRUE)
write.table(df_treatment, file='data_timeline_treatment_CMG_MM_formatted.txt', sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")

