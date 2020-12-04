library(purrr)
library(dplyr)
library(stringr)

setwd("C:/Users/abhmalat/OneDrive - Indiana University/cbio_MMRF")

df <- read.csv('raw_files/MMRF_IA15_summary.csv', sep=",", header=FALSE)

names(df) <- c("SAMPLE_ID","Patient","TISSUE","TIMEPOINT","CANCER_TYPE","t(4;14) cons",
               "t(11;14) cons","t(6;14) cons","t(14;16) cons","t(14;20) cons","t(8;14) cons",
               "None","HYPERDIPLOID","MYC","MUTATIONAL_SIGNATURE")

head(df,10)

df <- df[2:nrow(df), ]

df[df == '#N/A'] <- NA
df$TIMEPOINT[df$TIMEPOINT == 'ND'] <- 'Newly Diagnosed'
df$TISSUE[df$TISSUE == 'BM'] <- 'Bone Marrow'
df$TISSUE[df$TISSUE == 'PB'] <- 'Peripheral Blood'

df$TRANSLOCATION <- names(df[, c("t(4;14) cons","t(11;14) cons","t(6;14) cons",
                                 "t(14;16) cons","t(14;20) cons","t(8;14) cons")])[max.col(df[, c("t(4;14) cons", "t(11;14) cons","t(6;14) cons",
                                                                                                  "t(14;16) cons","t(14;20) cons","t(8;14) cons")],
                                                                                           ties.method = "first")]
df$TRANSLOCATION <- ifelse(df$None == 'None', 'None', df$TRANSLOCATION)

# Get patient ID from specimen ID
df$PATIENT_ID <- substr(df$SAMPLE_ID, 0, str_length(df$SAMPLE_ID)-5)


df$CANCER_TYPE_DETAILED <- ifelse(df$CANCER_TYPE == 'MM', 'Multiple Myeloma',
                                  ifelse(df$CANCER_TYPE == 'PCL', 'Plasma Cell Leukemia',
                                         df$CANCER_TYPE))

#df <- dplyr::select(df, c('PATIENT_ID', 'SAMPLE_ID', 'TISSUE', 'CANCER_TYPE',
#                          'CANCER_TYPE_DETAILED', 'TIMEPOINT', 'TRANSLOCATION',
#                          'HYPERDIPLOID','MYC','MUTATIONAL_SIGNATURE'))

df_old <- read.csv('data_clinical_sample.txt', sep="\t", header=FALSE, comment.char = "%")
colnames(df_old) <- df_old[5, ]

df_final <- df_old[1:5, ]
df_final <- rbind(df_final, df %>% select(c('PATIENT_ID', 'SAMPLE_ID', 'TISSUE', 'CANCER_TYPE',
                                      'CANCER_TYPE_DETAILED', 'TIMEPOINT', 'TRANSLOCATION')),
                                      df_old[6:nrow(df_old),])
# Remove duplicate samples
df_final <- unique(df_final, by = c('SAMPLE_ID', 'PATIENT_ID'))
# Fill blank values with NA
df_final[df_final == ''] <- NA

df_final <- left_join(df_final, df[, c('SAMPLE_ID', 'PATIENT_ID', 'HYPERDIPLOID','MYC','MUTATIONAL_SIGNATURE')],
                      by = c('SAMPLE_ID', 'PATIENT_ID'))
df_final[1, c('HYPERDIPLOID','MYC','MUTATIONAL_SIGNATURE')] <- c('Hyperdiploid','MYC','Mutational Signature')
df_final[2, c('HYPERDIPLOID','MYC','MUTATIONAL_SIGNATURE')] <- c('Hyperdiploid','MYC','Mutational Signature')
df_final[3, c('HYPERDIPLOID','MYC','MUTATIONAL_SIGNATURE')] <- c('STRING', 'STRING', 'STRING')
df_final[4, c('HYPERDIPLOID','MYC','MUTATIONAL_SIGNATURE')] <- c(8,9,10)
df_final[5, c('HYPERDIPLOID','MYC','MUTATIONAL_SIGNATURE')] <- c('HYPERDIPLOID','MYC','MUTATIONAL_SIGNATURE')

# ====================
# df_final <- data.table(c())
#df_patient_samples <- df %>% dplyr::select(PUBLIC_ID, SAMPLE_ID)
#
#names(df_patient_samples) <- c("PATIENT_ID", "SAMPLE_ID")
#
#summary(df_patient_samples)
#
#head(df_patient_samples)
#
outputFile <- "data_clinical_sample_1203.txt"

#file.remove(outputFile)
#file.create(outputFile)
#f <- file(outputFile)
#writeLines(c("#Patient Identifier\tSample Identifier ",
#             "#Patient Identifier\tSample Identifier",
#             "#STRING\tSTRING",
#             "#1\t2",
#             "PATIENT_ID\tSAMPLE_ID"), f)
#close(f)


write.table(df_final, file=outputFile, sep = "\t", col.names = FALSE,
            row.names=FALSE, quote = FALSE, append = FALSE)

metaFile <- file("meta_clinical_sample.txt")
writeLines(c(
  "cancer_study_identifier: mmrf_2020",
  "genetic_alteration_type: CLINICAL",
  "datatype: SAMPLE_ATTRIBUTES",
  paste("data_filename: ", outputFile)
), metaFile
)
close(metaFile)

print('Samples completed')
