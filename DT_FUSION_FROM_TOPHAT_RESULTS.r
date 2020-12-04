library(purrr)
library(dplyr)
library(org.Hs.eg.db)
library(annotate)
library(biomaRt)
library(tidyr)

setwd('c:/Users/abhmalat/OneDrive - Indiana University/MMRF_CoMMpass_IA16a')

df <- read.csv('fusion_transcripts/MMRF_CoMMpass_IA16a_TophatFusion_Results.txt', sep="\t", header=TRUE)

head(df)

colnames(df)

ensembl <- useEnsembl(host="grch37.ensembl.org",
                   biomart='ENSEMBL_MART_ENSEMBL',
                   dataset='hsapiens_gene_ensembl')


filters <- listFilters(mart=ensembl)

Hugo_Symbols <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id",
                      values=unique(c(df$left.Gene, df$right.Gene)), uniqueRows = TRUE, ensembl)
Entrez_Gene_Ids <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), filters = 'ensembl_gene_id',
                         values=unique(c(df$left.Gene, df$right.Gene)), uniqueRows = TRUE, ensembl)
Hugo_Symbols_name <- getBM(attributes = c("external_gene_name", "hgnc_symbol"), filters = "external_gene_name",
                      values=unique(c(df$left.Gene, df$right.Gene)), uniqueRows = TRUE, ensembl)
Entrez_Gene_Ids_name <- getBM(attributes = c("external_gene_name", "entrezgene_id"), filters = 'external_gene_name',
                         values=unique(c(df$left.Gene, df$right.Gene)), uniqueRows = TRUE, ensembl)


summary(Hugo_Symbols)
head(Hugo_Symbols)

summary(Entrez_Gene_Ids)
head(Entrez_Gene_Ids)

EntrezGrouped <- group_by(.data=Entrez_Gene_Ids, ensembl_gene_id)
EntrezGroupedName <- group_by(.data=Entrez_Gene_Ids_name, external_gene_name)
HugoGrouped <- group_by(.data=Hugo_Symbols, ensembl_gene_id)
HugoGroupedName <- group_by(.data=Hugo_Symbols_name, external_gene_name)

dfIntermediate <- df[, c('ID', 'left.Gene', 'right.Gene')] %>%
  left_join(EntrezGrouped, by=c("left.Gene" = "ensembl_gene_id")) %>%
  left_join(EntrezGrouped, by=c("right.Gene" = "ensembl_gene_id"))

dfIntermediate <- dfIntermediate %>%
  left_join(EntrezGroupedName, by=c("left.Gene" = "external_gene_name")) %>%
  left_join(EntrezGroupedName, by=c("right.Gene" = "external_gene_name"))

dfIntermediate <- dfIntermediate %>%
  left_join(HugoGrouped, by=c("left.Gene" = "ensembl_gene_id")) %>%
  left_join(HugoGrouped, by=c("right.Gene" = "ensembl_gene_id"))

dfIntermediate <- dfIntermediate %>%
  left_join(HugoGroupedName, by=c("left.Gene" = "external_gene_name")) %>%
  left_join(HugoGroupedName, by=c("right.Gene" = "external_gene_name"))

#dfIntermediate <- head(dfIntermediate, 100)
# _Left and _Right columns for Hugo and Entrez below are only used as intermediaries to create the Fusion label
dfIntermediate$Entrez_Left <- ifelse(is.na(dfIntermediate$entrezgene_id.x),
                                     dfIntermediate$left.Gene, dfIntermediate$entrezgene_id.x)
dfIntermediate$Entrez_Right <- ifelse(is.na(dfIntermediate$entrezgene_id.y),
                                      dfIntermediate$right.Gene, dfIntermediate$entrezgene_id.y)

dfIntermediate$Hugo_Left <- ifelse(is.na(dfIntermediate$hgnc_symbol.x),
                                   dfIntermediate$left.Gene, dfIntermediate$hgnc_symbol.x)
dfIntermediate$Hugo_Right <- ifelse(is.na(dfIntermediate$hgnc_symbol.y),
                                    dfIntermediate$right.Gene, dfIntermediate$hgnc_symbol.y)

dfIntermediate$Fusion <- paste(dfIntermediate$Hugo_Left, '-', dfIntermediate$Hugo_Right, ' Fusion')

dfIntermediate$center <- 'sequencing_center_unknown'
dfIntermediate$DNA_support <- 'no'
dfIntermediate$RNA_support <- 'yes'
dfIntermediate$Method <- 'Tophat-fusion'
dfIntermediate$Frame <- 'Unknown'

dfLeft <- dfIntermediate[,c('left.Gene','hgnc_symbol.x',  'entrezgene_id.x', 'center', 'ID', 'Fusion', 'DNA_support', 'RNA_support', 'Method', 'Frame')]
dfRight <- dfIntermediate[,c('right.Gene', 'hgnc_symbol.y',  'entrezgene_id.y', 'center', 'ID', 'Fusion', 'DNA_support', 'RNA_support', 'Method', 'Frame')]

names(dfLeft) <- c('TophatResultsGeneName', 'Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'Tumor_Sample_Barcode', 'Fusion', 'DNA_support', 'RNA_support', 'Method', 'Frame')
names(dfRight) <- c('TophatResultsGeneName', 'Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'Tumor_Sample_Barcode', 'Fusion', 'DNA_support', 'RNA_support', 'Method', 'Frame')

dfFinal <- rbind(dfLeft, dfRight)

outputFile <- "data_fusion_grch37.txt"
# outputFileWithBlanks <- "zdata_fusion_blanks.txt"

if(file.exists(outputFile)) {
    file.remove(outputFile)
}
file.create(outputFile)

write.table(dfFinal %>%
              filter(dfFinal$Hugo_Symbol != "" & dfFinal$Entrez_Gene_Id != "") %>%
              select(c('Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'Tumor_Sample_Barcode', 'Fusion', 'DNA_support', 'RNA_support', 'Method', 'Frame')),
            file=outputFile, sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, append = TRUE, na="NA")

# write.table(dfFinal, file=outputFileWithBlanks, sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, append = TRUE, na="NA")

print('Completed.')