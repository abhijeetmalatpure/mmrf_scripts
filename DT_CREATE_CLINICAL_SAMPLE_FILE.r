library(purrr)
library(dplyr)

setwd("f:/mmrf_data")

df <- read.csv('README_files/MMRF_CoMMpass_IA11_Summary_Table.txt', sep="\t", header=FALSE)

names(df) <- c("PUBLIC_ID", "SPECTRUM_SEQ", "SAMPLE_ID", "ANALYSIS_TYPE")

head(df,10)

summary(df)

df_patient_samples <- df %>% dplyr::select(PUBLIC_ID, SAMPLE_ID)

names(df_patient_samples) <- c("PATIENT_ID", "SAMPLE_ID")

summary(df_patient_samples)

head(df_patient_samples)

outputFile <- "DT_R/cbio_files/data_clinical_sample.txt"

file.remove(outputFile)
file.create(outputFile)
f <- file(outputFile)
writeLines(c("#Patient Identifier\tSample Identifier ",
             "#Patient Identifier\tSample Identifier",
             "#STRING\tSTRING",
             "#1\t2",
             "PATIENT_ID\tSAMPLE_ID"), f)
close(f)


write.table(distinct(df_patient_samples[,c("PATIENT_ID","SAMPLE_ID")]), file=outputFile, sep = "\t", col.names = FALSE, row.names=FALSE, quote = FALSE, append = TRUE) 

