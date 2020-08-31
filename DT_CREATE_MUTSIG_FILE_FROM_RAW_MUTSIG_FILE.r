library(purrr)
library(dplyr)

setwd("f:/mmrf_data")

df <- read.csv('DT_R/cbio_files/MMRF.sig_genes.txt', sep="\t", header=TRUE)

#names(df) = c("PUBLIC_ID", "SPECTRUM_SEQ", "SAMPLE_ID", "ANALYSIS_TYPE")

head(df,10)

summary(df)

df$rank <- seq.int(nrow(df))

df_cbio_mutsig <- df %>% dplyr::select(rank, gene, N_nonsilent, n_nonsilent, p, q)

names(df_cbio_mutsig) <- c("rank", "gene", "N", "n", "p", "q")

summary(df_cbio_mutsig)

head(df_cbio_mutsig, 100)

outputFile <- "DT_R/cbio_files/data_mutsig.txt"

file.remove(outputFile)

write.table(df_cbio_mutsig, file=outputFile, sep = "\t", col.names = FALSE, row.names=FALSE, quote = FALSE, append = FALSE) 

