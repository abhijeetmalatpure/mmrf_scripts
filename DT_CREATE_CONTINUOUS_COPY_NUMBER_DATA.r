library(purrr)
library(org.Hs.eg.db)
library(annotate)
library(dplyr)
library(tidyr)

library(org.Hs.eg.db)
library(annotate)
library(biomaRt)

setwd('c:/Users/abhmalat/OneDrive - Indiana University/MMRF_CoMMpass_IA16a')

df <- read.csv('copy_number_estimates/MMRF_CoMMpass_IA16a_CNA_LongInsert_PerGene_LowestSegment.txt', sep="\t", header=TRUE)

head(df)
colnames(df)
summary(df)

# Lookup Gene by Gene (ensembl id)
ensembl <- useMart(host='apr2020.archive.ensembl.org',
                   biomart='ENSEMBL_MART_ENSEMBL',
                   dataset='hsapiens_gene_ensembl')

Hugo_Symbols <- getBM(c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values=df$Gene, ensembl)

Entrez_Gene_Ids <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), filters = 'ensembl_gene_id', values=df$Gene, ensembl)

Entrez_Gene_Ids[,"ensembl_gene_id"] <- as.factor(Entrez_Gene_Ids[, "ensembl_gene_id"])

Entrez_Gene_Ids[,"entrezgene_id"] <- as.factor(Entrez_Gene_Ids[, "entrezgene_id"])

foo <- unique(arrange(Entrez_Gene_Ids, ensembl_gene_id, entrezgene_id))

fooHugo <- left_join(foo, Hugo_Symbols, by="ensembl_gene_id")

head(fooHugo,100)

nrow(fooHugo)

summary(fooHugo)

dfEntrez <- left_join(df, fooHugo, by=c("Gene" = "ensembl_gene_id") )

nrow(dfEntrez)

colnames(dfEntrez)

levels ( dfEntrez$entrezgene_id )

dfWithHugo <- dfEntrez %>% dplyr::select(hgnc_symbol, entrezgene_id, everything())

colnames(dfWithHugo[,1:5])

names(dfWithHugo)[1] <- "Hugo_Symbol"
names(dfWithHugo)[2] <- "Entrez_Gene_Id"

colnames(dfWithHugo)

dfPrelim <- dfWithHugo %>%
  filter(Hugo_Symbol != "") %>%
  dplyr::select(-c(Gene))

head(dfPrelim[,1:8],300)

outputFile <- ("data_CNA_LongInsert_PerGene_Lowest_Continuous.txt")

if(file.exists(outputFile)) {
    file.remove(outputFile)
}
file.create(outputFile)

write.table(dfPrelim, file=outputFile, sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, append = TRUE, na="NA")

# output case_list for _cna
sampleIds <- colnames(dfPrelim[,3:ncol(dfPrelim)])

#outputFileCL <- ("case_lists/cases_cna_lowest_continuous.txt")
#
#if(file.exists(outputFileCL)) {
#    file.remove(outputFileCL)
#}
#file.create(outputFileCL)

#f <- file(outputFileCL)
#writeLines(c(
#  "cancer_study_identifier: mmrf_2020",
#  "stable_id: mmrf_2020_cna_lowest",
#  "case_list_name: Continuous Copy Number Analysis Lowest Segment",
#  "case_list_description: Continuous Copy Number Analysis Lowest Segment",
#  paste("case_list_ids: ", paste(sampleIds, collapse = '\t'))
#), f
#)
#close(f)
print("Completed.")
