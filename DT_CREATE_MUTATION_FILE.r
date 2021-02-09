library(purrr)
library(org.Hs.eg.db)
library(annotate)
library(dplyr)
library(tidyr)
library(biomaRt)

setwd('c:/Users/abhmalat/OneDrive - Indiana University/cbio_MMRF')

df <- read.csv('somatic_mutation_files/MMRF_CoMMpass_IA16a_All_Canonical_NS_Variants.txt', sep="\t", header=TRUE)

head(df)
colnames(df)
summary(df)

count(df,ANN....EFFECT, sort = TRUE)


# SEPARATE ANN....EFFECT field into individual columns
df <- separate(data=df, col=ANN....EFFECT, into = c("ANN_EFFECT1","ANN_EFFECT2", "ANN_EFFECT3", "ANN_EFFECT4"), sep="&", remove=FALSE, fill= "right")

levels(df$ANN....BIOTYPE)

ensembl <- useMart(host="apr2020.archive.ensembl.org",
                   biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')

Entrez_Gene_Ids <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), filters = 'ensembl_gene_id', values=df$ANN....GENEID, ensembl)

Entrez_Gene_Ids[,"ensembl_gene_id"] <- as.factor(Entrez_Gene_Ids[,"ensembl_gene_id"])

Entrez_Gene_Ids[,"entrezgene_id"] <- as.factor(Entrez_Gene_Ids[,"entrezgene_id"])

foo <- unique(arrange(Entrez_Gene_Ids, ensembl_gene_id, entrezgene_id))

head(foo,100)

nrow(foo)

summary(Entrez_Gene_Ids)

class(Entrez_Gene_Ids[,"ensembl_gene_id"])

levels(Entrez_Gene_Ids[,"ensembl_gene_id"])

levels(foo[,"ensembl_gene_id"])

nrow(Entrez_Gene_Ids)

class(df[,"ANN....GENEID"])

head(Entrez_Gene_Ids, 10)

nrow(df)

levels(df[,"ANN....GENEID"])

dfEntrez <- left_join(df, foo, by=c("ANN....GENEID" = "ensembl_gene_id") )

nrow(dfEntrez)

colnames(dfEntrez)

dfVariants <- data.frame("ANN_EFFECT1" = c("disruptive_inframe_insertion", "inframe_insertion", "disruptive_inframe_deletion", "inframe_deletion", "missense_variant", "stop_gained", "synonymous_variant", "stop_retained_variant", "splice_donor_variant", "splice_acceptor_variant", "start_lost", "initiator_codon_variant", "stop_lost", "3_prime_UTR_variant", "downstream_gene_variant", "5_prime_UTR_variant", "5_prime_UTR_premature_start_codon_gain_variant", "5_prime_UTR_truncation", "upstream_gene_variant", "non_coding_exon_variant", "exon_loss_variant", "intron_variant", "intragenic_variant", "splice_region_variant", "frameshift_variant"),
                         "Variant_Classification" = c("In_Frame_Ins", "In_Frame_Ins", "In_Frame_Del", "In_Frame_Del", "Missense_Mutation", "Nonsense_Mutation", "Silent", "Silent", "Splice_Site", "Splice_Site", "Translation_Start_Site", "Translation_Start_Site", "Nonstop_Mutation", "3'UTR", "3'Flank", "5'UTR", "5'UTR", "5'UTR", "5'Flank", "Intron", "Splice_Site", "Intron", "Intron", "Splice_Region", "Unknown"))
head(dfVariants, 50)

dfVarEntrez <- left_join(dfEntrez, dfVariants, by="ANN_EFFECT1")

head(dfVarEntrez)

dfVarEntrez2 <- dfVarEntrez %>% dplyr::select( ANN....GENE, entrezgene_id, X.CHROM, POS, Variant_Classification, REF, ALT, ID, Sample, ANN....HGVS_P, ANN....AA_POS )

head(dfVarEntrez2)
names(dfVarEntrez2) <- c("Hugo_Symbol", "Entrez_Gene_Id", "Chromosome", "Start_Position", "Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele1", "dbSNP_RS", "Tumor_Sample_Barcode", "HGVSp_Short", "Protein_Pos")

head(dfVarEntrez2,30)

aaa0<-as.matrix(dfVarEntrez2)

ddd<-nchar(aaa0[,"Reference_Allele"])-nchar(aaa0[,"Tumor_Seq_Allele1"])
ddd[which(ddd<0)]<-0
ddd1<-nchar(aaa0[,"Reference_Allele"])-nchar(aaa0[,"Tumor_Seq_Allele1"])
Start_Position<-as.numeric(aaa0[,"Start_Position"])
End_Position<-as.numeric(aaa0[,"Start_Position"])+ddd

Start_Position
End_Position

eee<-rep("",length(ddd1))
eee[which(ddd1==0)]<-"SNP"
eee[which(ddd1>0)]<-"DEL"
eee[which(ddd1<0)]<-"INS"

#eee

dfVarEntrez2$End_Position = End_Position

dfVarEntrez2$Variant_Type = eee

dfVarEntrez2$NCBI_Build='GRCh37'

colnames(dfVarEntrez2)

output_cols <- c(
  "Hugo_Symbol", "Entrez_Gene_Id", "Chromosome",
  "Start_Position", "End_Position",
  "Variant_Classification", "Variant_Type",
  "Reference_Allele", "Tumor_Seq_Allele1",
  "dbSNP_RS", "Tumor_Sample_Barcode",
  "HGVSp_Short", "Protein_Pos"
)

dfVarEntrez3 <- dfVarEntrez2 %>% dplyr::select("Hugo_Symbol", "Entrez_Gene_Id", "Chromosome",
                                               "Start_Position", "End_Position",
                                               "Variant_Classification", "Variant_Type",
                                               "Reference_Allele", "Tumor_Seq_Allele1",
                                               "dbSNP_RS", "Tumor_Sample_Barcode",
                                               "HGVSp_Short", "Protein_Pos")

outputFile <- ("somatic_mutation_files/data_mutation_file_37.maf")
if(file.exists(outputFile)) {
    file.remove(outputFile)
}
file.create(outputFile)

write.table(dfVarEntrez3, file=outputFile, sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, append = TRUE, na="")

# output case_list for _sequence
sampleIds <- str(dfVarEntrez2$Tumor_Sample_Barcode, )

outputFileCL <- ("case_lists/cases_sequenced.txt")


if(file.exists(outputFileCL)) {
    file.remove(outputFileCL)
}
file.create(outputFileCL)

f <- file(outputFileCL)
writeLines(c(
  "cancer_study_identifier: mmrf_2020",
  "stable_id: mmrf_2020_sequenced",
  "case_list_name: Sequenced Tumors",
  "case_list_description: All Sequenced Tumors",
  paste("case_list_ids: ", paste(dfVarEntrez2$Tumor_Sample_Barcode, collapse = '\t'))
  ), f
)
print('Completed.')
close(f)