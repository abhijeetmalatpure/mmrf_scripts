library(purrr)
library(dplyr)
library(org.Hs.eg.db)
library(annotate)

library(biomaRt)

setwd('c:/Users/abhmalat/OneDrive - Indiana University/MMRF_CoMMpass_IA16a')

df <- read.csv('gene_expression_estimates/MMRF_CoMMpass_IA16a_E74GTF_Sailfish_Gene_TPM.txt', sep="\t", header=TRUE)

head(df)

colnames(df)

ensembl <- useEnsembl(host="grch37.ensembl.org",
                   biomart='ENSEMBL_MART_ENSEMBL',
                   dataset='hsapiens_gene_ensembl')


listFilters(mart=ensembl)

Hugo_Symbols <- getBM(c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values=df$GENE_ID, ensembl)

Entrez_Gene_Ids <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), filters = 'ensembl_gene_id', values=df$GENE_ID, ensembl)


summary(Hugo_Symbols)

head(Hugo_Symbols)

summary(Entrez_Gene_Ids)

head(Entrez_Gene_Ids)

EntrezGrouped <- group_by(.data=Entrez_Gene_Ids, ensembl_gene_id)

dfEntrez <- left_join(df, EntrezGrouped, by=c("GENE_ID" = "ensembl_gene_id"))

dfEntrez <- dfEntrez %>% dplyr::select( entrezgene_id, everything() )

head(dfEntrez[,1:5])

dfHugo <- left_join(dfEntrez, Hugo_Symbols, by=c("GENE_ID" = "ensembl_gene_id"))

dfHugo <- dfHugo %>% dplyr::select(hgnc_symbol, everything())

head(dfHugo[,1:6])

dfFinal <- dfHugo %>% 
            filter(hgnc_symbol != "") %>%
            dplyr::select(-GENE_ID)

nrow(dfFinal)


head(dfFinal[,1:6])

names(dfFinal)[1] <- "Hugo_Symbol"

names(dfFinal)[2] <- "Entrez_Gene_Id"

head(dfFinal[,1:6])  


outputFile <- "data_expression_rnaseq_tpm_grch37.txt"

if(file.exists(outputFile)) {
    file.remove(outputFile)
}
file.create(outputFile)

write.table(dfFinal, file=outputFile, sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, append = TRUE, na="NA")


# output case_list for _

sampleIds <- colnames(dfFinal[,3:ncol(dfFinal)])

outputFileCL <- "case_lists/cases_rna_seq_mrna.txt"

if(file.exists(outputFileCL)) {
    file.remove(outputFileCL)
}
file.create(outputFileCL)

f <- file(outputFileCL)
writeLines(
  c(
  "cancer_study_identifier: mmrf_2020",
  "stable_id: mmrf_2020_rna_seq_mrna",
  "case_list_name: RNA Seq",
  "case_list_description: E74GTF_Sailfish_Gene_TPM",
  paste("case_list_ids: ", paste(sampleIds, collapse = '\t'))
  ), f)
close(f)
print('Completed.')