library(purrr)
library(org.Hs.eg.db)
library(annotate)
library(dplyr)
library(tidyr)

library(org.Hs.eg.db)
library(annotate)

setwd('f:/mmrf_data')

df <- read.csv('copy_number_files/MMRF_CoMMpass_IA11a_CNA_Exome_extended_PerGene_LargestSegment.txt', sep="\t", header=TRUE)

head(df)
colnames(df)
summary(df)

# Lookup Gene by Gene (ensembl id)
ensembl <- useMart(host='dec2016.archive.ensembl.org',
                   biomart='ENSEMBL_MART_ENSEMBL',
                   dataset='hsapiens_gene_ensembl')

Hugo_Symbols <- getBM(c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values=df$Gene, ensembl)

Entrez_Gene_Ids <- getBM(attributes = c("ensembl_gene_id", "entrezgene"), filters = 'ensembl_gene_id', values=df$Gene, ensembl)

Entrez_Gene_Ids[,"ensembl_gene_id"] <- as.factor(Entrez_Gene_Ids[, "ensembl_gene_id"])

Entrez_Gene_Ids[,"entrezgene"] <- as.factor(Entrez_Gene_Ids[, "entrezgene"])

foo <- unique(arrange(Entrez_Gene_Ids, ensembl_gene_id, entrezgene))

fooHugo <- left_join(foo, Hugo_Symbols, by="ensembl_gene_id")

head(fooHugo,100)

nrow(fooHugo)

summary(fooHugo)

dfEntrez <- left_join(df, fooHugo, by=c("Gene" = "ensembl_gene_id") )

nrow(dfEntrez)

colnames(dfEntrez)

levels ( dfEntrez$entrezgene )

dfWithHugo <- dfEntrez %>% dplyr::select(hgnc_symbol, entrezgene, everything())

colnames(dfWithHugo[,1:5])

names(dfWithHugo)[1] <- "Hugo_Symbol"
names(dfWithHugo)[2] <- "Entrez_Gene_Id"

colnames(dfWithHugo)

dfPrelim <- dfWithHugo %>% 
  filter(Hugo_Symbol != "") %>%
  dplyr::select(-Gene)


head(dfPrelim[,1:10])

calcCopyNumLevel <- function(x) {
  if ( x < -2 ) {
    return(-2)
  } else
  if ( x < -0.3) {
    return(-1)
  } else
  if ( x < 0.3 ) {
    return(0)
  } else
  if ( x < 2 ) {
    return(1)
  } else {
    return(2)
  }
}

dfPrelim[,3:ncol(dfPrelim)] = as.data.frame(lapply(dfPrelim[,3:ncol(dfPrelim)], FUN= function(x) sapply(x, FUN=calcCopyNumLevel) ) )

head(dfPrelim[,1:8],300)

outputFile <- "DT_R/cbio_files/data_CNA.txt"

file.remove(outputFile)
file.create(outputFile)

write.table(dfPrelim, file=outputFile, sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, append = TRUE, na="") 

# output case_list for _cna
sampleIds <- colnames(dfPrelim[,3:ncol(dfPrelim)])

outputFileCL <- "DT_R/cbio_files/case_lists/cases_cna.txt"

file.remove(outputFileCL)
file.create(outputFileCL)

f <- file(outputFileCL)
writeLines(c(
  "cancer_study_identifier: mmrf_mar_2018",
  "stable_id: mmrf_mar_2018_cna",
  "case_list_name: Discrete Copy Number Analysis",
  "case_list_description: Discrete Copy Number Analysis",
  paste("case_list_ids: ", paste(sampleIds, collapse = '\t'))
), f
)
close(f)


