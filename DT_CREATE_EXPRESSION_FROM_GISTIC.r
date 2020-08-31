library(purrr)
library(dplyr)
library(org.Hs.eg.db)
library(annotate)

setwd('f:/mmrf_data')

df <- read.csv('expression_estimates_gene-based/MMRF_CoMMpass_IA11a_E74cDNA_Salmon_Gene_TPM.txt', sep="\t", header=TRUE)

head(df)

colnames(df)

ensembl <- useMart(host='dec2016.archive.ensembl.org',
                   biomart='ENSEMBL_MART_ENSEMBL',
                   dataset='hsapiens_gene_ensembl')


listFilters(mart=ensembl)

Hugo_Symbols <- getBM(c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values=df$GENE_ID, ensembl)

Entrez_Gene_Ids <- getBM(attributes = c("ensembl_gene_id", "entrezgene"), filters = 'ensembl_gene_id', values=df$GENE_ID, ensembl)


summary(Hugo_Symbols)

head(Hugo_Symbols)

summary(Entrez_Gene_Ids)

head(Entrez_Gene_Ids)

EntrezGrouped <- group_by(.data=Entrez_Gene_Ids, ensembl_gene_id)

dfEntrez <- left_join(df, EntrezGrouped, by=c("GENE_ID" = "ensembl_gene_id"))

colnames(dfEntrez)[866]

dfEntrez <- dfEntrez %>% dplyr::select( entrezgene, everything() )

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


outputFile <- "DT_R/cbio_files/data_expression_file.txt"

file.remove(outputFile)
file.create(outputFile)

write.table(dfFinal, file=outputFile, sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, append = TRUE, na="") 


# output case_list for _

sampleIds <- colnames(dfFinal)

outputFileCL <- "DT_R/cbio_files/case_lists/cases_rna_seq_mrna.txt"

file.remove(outputFileCL)
file.create(outputFileCL)

f <- file(outputFileCL)
writeLines(
  c(
  "cancer_study_identifier: mmrf_mar_2018",
  "stable_id: mmrf_mar_2018_rna_seq_mrna",
  "case_list_name: RNA Seq",
  "case_list_description: All Expression Levels TPN",
  paste("case_list_ids: ", paste(sampleIds[3:866], collapse = '\t'))
  ), f)
close(f)
