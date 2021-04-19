ncol(df[29252,])
df <- read.csv('data_timeline_lab_test_CMG_MM.txt', sep="\t", header=TRUE, quote = "escape")

df <- read.csv('data_timeline_lab_test_CMG_MM.txt', sep="\t", header=TRUE, quote = "")

df[29251, ]
write.table(df, file='data_timeline_lab_test_CMG_MM_formatted.txt', sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")

library(purrr)
library(dplyr)
library(org.Hs.eg.db)
library(annotate)

setwd('C:/Users/abhmalat/OneDrive - Indiana University/cbio_MM/mm_files_1122b')


# Fix quotes for height
df <- read.csv('data_timeline_lab_test_CMG_MM.txt', sep="\t", header=TRUE, quote = "")

setwd('C:/Users/abhmalat/OneDrive - Indiana University/cbio_MM/')

# Fix quotes for height
df <- read.csv('data_timeline_lab_test_CMG_MM.txt', sep="\t", header=TRUE, quote = "")
unique(df$grade)
unique(df$total_cells)
unique(df$nucleic_acid_derivative)
unique(df$parent_specimen_availability)
unique(df$specimen_status)


# Fix quotes for height
df <- read.csv('data_timeline_lab_test_CMG_MM.txt', sep="\t", header=TRUE, quote = "")
write.table(df, file='data_timeline_lab_test_CMG_MM_formatted.txt', sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")


df_treatment <- read.csv('data_timeline_treatment_CMG_MM.txt', sep="\t", header=TRUE)
write.table(df_treatment, file='data_timeline_treatment_CMG_MM_formatted.txt', sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")

library(purrr)
library(dplyr)

setwd("C:\Users\abhmalat\OneDrive - Indiana University\cbio_MMRF\raw_files")

setwd("C:/Users/abhmalat/OneDrive - Indiana University/cbio_MMRF/raw_files")

df <- read.csv('MMRF_IA15_summary.csv', sep="\t", header=TRUE)

df <- read.csv('MMRF_IA15_summary.csv', sep=",", header=TRUE)

df <- read.csv('MMRF_IA15_summary.csv', sep=",", header=FALSE)


names(df) <- c("Specimen_ID","Patient","Tissue","Timepoint","Disease","t(4;14) cons",
               "t(11;14) cons","t(6;14) cons","t(14;16) cons","t(14;20) cons","t(8;14) cons",
               "None","Hyperdiploid","MYC","Mutational Signature")

head(df,10)

for (col in colnames(df)) set(df, i=which(df[[col]] =='#N/A'), j=col, value="NA")

setwd("C:/Users/abhmalat/OneDrive - Indiana University/cbio_MMRF")

df <- read.csv('raw_files/MMRF_IA15_summary.csv', sep=",", header=FALSE)


names(df) <- c("Specimen_ID","Patient","Tissue","Timepoint","Disease","t(4;14) cons",
               "t(11;14) cons","t(6;14) cons","t(14;16) cons","t(14;20) cons","t(8;14) cons",
               "None","Hyperdiploid","MYC","Mutational Signature")

head(df,10)

df <- df[2:nrow(df), ]

library(stringr)

df$PATIENT_ID <- str_extract(df$Specimen_ID, '$[_[0-9]_')

df$PATIENT_ID <- str_extract(df$Specimen_ID, '$_[0-9]_')

df$PATIENT_ID <- str_extract(df$Specimen_ID, '$[_][0-9][_]')

df$PATIENT_ID <- str_extract(df$Specimen_ID, '${_}[0-9]{_}')

df$PATIENT_ID <- str_extract(df$Specimen_ID, '^.[_][0-9][_]')

df$PATIENT_ID <- substr(df$Specimen_ID, 0, str_length(df$Specimen_ID)-5)
library(purrr)
library(dplyr)
library(stringr)

setwd("C:/Users/abhmalat/OneDrive - Indiana University/cbio_MMRF")

df <- read.csv('raw_files/MMRF_IA15_summary.csv', sep=",", header=FALSE)

names(df) <- c("SAMPLE_ID","Patient","Tissue","Timepoint","Disease","t(4;14) cons",
               "t(11;14) cons","t(6;14) cons","t(14;16) cons","t(14;20) cons","t(8;14) cons",
               "None","Hyperdiploid","MYC","Mutational Signature")

head(df,10)

df <- df[2:nrow(df), ]

df[df == '#N/A'] <- 'NA'

# Get patient ID from specimen ID
df$PATIENT_ID <- substr(df$Specimen_ID, 0, str_length(df$Specimen_ID)-5)

df$PATIENT_ID <- substr(df$SAMPLE_ID, 0, str_length(df$SAMPLE_ID)-5)

library(purrr)
library(dplyr)
library(stringr)

setwd("C:/Users/abhmalat/OneDrive - Indiana University/cbio_MMRF")

df <- read.csv('raw_files/MMRF_IA15_summary.csv', sep=",", header=FALSE)

names(df) <- c("SAMPLE_ID","Patient","Tissue","Timepoint","Disease","t(4;14) cons",
               "t(11;14) cons","t(6;14) cons","t(14;16) cons","t(14;20) cons","t(8;14) cons",
               "None","Hyperdiploid","MYC","Mutational Signature")

head(df,10)

df <- df[2:nrow(df), ]

df[df == '#N/A'] <- 'NA'

# Get patient ID from specimen ID
df$PATIENT_ID <- substr(df$SAMPLE_ID, 0, str_length(df$SAMPLE_ID)-5)

summary(df)
df <- read.csv('data_clinical_sample.txt', sep="\t", header=TRUE)


df <- read.csv('raw_files/MMRF_IA15_summary.csv', sep=",", header=FALSE)

names(df) <- c("SAMPLE_ID","Patient","Tissue","Timepoint","Disease","t(4;14) cons",
               "t(11;14) cons","t(6;14) cons","t(14;16) cons","t(14;20) cons","t(8;14) cons",
               "None","Hyperdiploid","MYC","Mutational Signature")

head(df,10)

df <- df[2:nrow(df), ]

df[df == '#N/A'] <- 'NA'

# Get patient ID from specimen ID
df$PATIENT_ID <- substr(df$SAMPLE_ID, 0, str_length(df$SAMPLE_ID)-5)

summary(df)
toupper(df[max.col(data), c("t(4;14) cons",
               "t(11;14) cons","t(6;14) cons","t(14;16) cons","t(14;20) cons","t(8;14) cons")])
toupper(names(df)[max.col(data), c("t(4;14) cons",
               "t(11;14) cons","t(6;14) cons","t(14;16) cons","t(14;20) cons","t(8;14) cons")])
toupper(names(df)[max.col(data)])


data$a <- c(1,0,0,1,0)
data$b <- c(0,1,1,0,0)
data$c <- c(0,0,0,0,1)

dt <- data.table()
data <- data.table()
data$a <- c(1,0,0,1,0)
data$b <- c(0,1,1,0,0)
data$c <- c(0,0,0,0,1)
data$tranformed <- toupper(names(data)[max.col(data)])
names(data)
names(data)[max.col(data)]
data <- data.table()
> data$a <- c(1,0,0,1,0)
> data$b <- c(0,1,1,0,0)
> data$c <- c(0,0,0,0,1)

data <- data.table('a' <- c(1,0,0,1,0), 'b' <- c(0,1,1,0,0), 'c' <-  c(0,0,0,0,1))
data
toupper(names(data))
data[max.col(data)]
max.col(data)
max.col(df[, c("t(4;14) cons",
               "t(11;14) cons","t(6;14) cons","t(14;16) cons","t(14;20) cons","t(8;14) cons")])
names(df[, c("t(4;14) cons",
               "t(11;14) cons","t(6;14) cons","t(14;16) cons","t(14;20) cons","t(8;14) cons")])[max.col(df[, c("t(4;14) cons",
               "t(11;14) cons","t(6;14) cons","t(14;16) cons","t(14;20) cons","t(8;14) cons")])]
df$Translocation <- names(df[, c("t(4;14) cons","t(11;14) cons","t(6;14) cons",
                                 "t(14;16) cons","t(14;20) cons","t(8;14) cons")])[max.col(df[, c("t(4;14) cons", "t(11;14) cons","t(6;14) cons",
                                                                                                  "t(14;16) cons","t(14;20) cons","t(8;14) cons")])]

library(purrr)
library(dplyr)
library(stringr)

setwd("C:/Users/abhmalat/OneDrive - Indiana University/cbio_MMRF")

df <- read.csv('raw_files/MMRF_IA15_summary.csv', sep=",", header=FALSE)

names(df) <- c("SAMPLE_ID","Patient","Tissue","Timepoint","Disease","t(4;14) cons",
               "t(11;14) cons","t(6;14) cons","t(14;16) cons","t(14;20) cons","t(8;14) cons",
               "None","Hyperdiploid","MYC","Mutational Signature")

head(df,10)

df <- df[2:nrow(df), ]

df[df == '#N/A'] <- 'NA'
df$Translocation <- names(df[, c("t(4;14) cons","t(11;14) cons","t(6;14) cons",
                                 "t(14;16) cons","t(14;20) cons","t(8;14) cons")])[max.col(df[, c("t(4;14) cons", "t(11;14) cons","t(6;14) cons",
                                                                                                  "t(14;16) cons","t(14;20) cons","t(8;14) cons")])]
df$Translocation <- which(df$None == 'None', 'None')
# Get patient ID from specimen ID
df$PATIENT_ID <- substr(df$SAMPLE_ID, 0, str_length(df$SAMPLE_ID)-5)

df_old <- read.csv('data_clinical_sample.txt', sep="\t", header=TRUE)
df$Translocation <- ifelse(df$None == 'None', 'None', df$Translocation)

library(purrr)
library(dplyr)
library(stringr)

setwd("C:/Users/abhmalat/OneDrive - Indiana University/cbio_MMRF")

df <- read.csv('raw_files/MMRF_IA15_summary.csv', sep=",", header=FALSE)

names(df) <- c("SAMPLE_ID","Patient","Tissue","Timepoint","Disease","t(4;14) cons",
               "t(11;14) cons","t(6;14) cons","t(14;16) cons","t(14;20) cons","t(8;14) cons",
               "None","Hyperdiploid","MYC","Mutational Signature")

head(df,10)

df <- df[2:nrow(df), ]

df[df == '#N/A'] <- 'NA'
df$Translocation <- names(df[, c("t(4;14) cons","t(11;14) cons","t(6;14) cons",
                                 "t(14;16) cons","t(14;20) cons","t(8;14) cons")])[max.col(df[, c("t(4;14) cons", "t(11;14) cons","t(6;14) cons",
                                                                                                  "t(14;16) cons","t(14;20) cons","t(8;14) cons")],
                                                                                           ties.method = "first")]
df$Translocation <- ifelse(df$None == 'None', 'None', df$Translocation)
# Get patient ID from specimen ID
df$PATIENT_ID <- substr(df$SAMPLE_ID, 0, str_length(df$SAMPLE_ID)-5)

df_old <- read.csv('data_clinical_sample.txt', sep="\t", header=TRUE)

df <- read.csv('raw_files/MMRF_IA15_summary.csv', sep=",", header=FALSE)

names(df) <- c("SAMPLE_ID","Patient","Tissue","Timepoint","CANCER_TYPE","t(4;14) cons",
               "t(11;14) cons","t(6;14) cons","t(14;16) cons","t(14;20) cons","t(8;14) cons",
               "None","Hyperdiploid","MYC","Mutational Signature")

library(purrr)
library(dplyr)
library(stringr)

setwd("C:/Users/abhmalat/OneDrive - Indiana University/cbio_MMRF")

df <- read.csv('raw_files/MMRF_IA15_summary.csv', sep=",", header=FALSE)

names(df) <- c("SAMPLE_ID","Patient","Tissue","Timepoint","CANCER_TYPE","t(4;14) cons",
               "t(11;14) cons","t(6;14) cons","t(14;16) cons","t(14;20) cons","t(8;14) cons",
               "None","Hyperdiploid","MYC","Mutational Signature")


df <- df[2:nrow(df), ]

df[df == '#N/A'] <- 'NA'

df$TRANSLOCATION <- names(df[, c("t(4;14) cons","t(11;14) cons","t(6;14) cons",
                                 "t(14;16) cons","t(14;20) cons","t(8;14) cons")])[max.col(df[, c("t(4;14) cons", "t(11;14) cons","t(6;14) cons",
                                                                                                  "t(14;16) cons","t(14;20) cons","t(8;14) cons")],
                                                                                           ties.method = "first")]
df$TRANSLOCATION <- ifelse(df$None == 'None', 'None', df$TRANSLOCATION)

df <- read.csv('raw_files/MMRF_IA15_summary.csv', sep=",", header=FALSE)

names(df) <- c("SAMPLE_ID","Patient","Tissue","TIMEPOINT","CANCER_TYPE","t(4;14) cons",
               "t(11;14) cons","t(6;14) cons","t(14;16) cons","t(14;20) cons","t(8;14) cons",
               "None","Hyperdiploid","MYC","Mutational Signature")

df <- read.csv('raw_files/MMRF_IA15_summary.csv', sep=",", header=FALSE)

names(df) <- c("SAMPLE_ID","Patient","TISSUE","TIMEPOINT","CANCER_TYPE","t(4;14) cons",
               "t(11;14) cons","t(6;14) cons","t(14;16) cons","t(14;20) cons","t(8;14) cons",
               "None","Hyperdiploid","MYC","Mutational Signature")

df <- df[2:nrow(df), ]

df[df == '#N/A'] <- 'NA'
df$TRANSLOCATION <- names(df[, c("t(4;14) cons","t(11;14) cons","t(6;14) cons",
                                 "t(14;16) cons","t(14;20) cons","t(8;14) cons")])[max.col(df[, c("t(4;14) cons", "t(11;14) cons","t(6;14) cons",
                                                                                                  "t(14;16) cons","t(14;20) cons","t(8;14) cons")],
                                                                                           ties.method = "first")]
df$TRANSLOCATION <- ifelse(df$None == 'None', 'None', df$TRANSLOCATION)

# Get patient ID from specimen ID
df$PATIENT_ID <- substr(df$SAMPLE_ID, 0, str_length(df$SAMPLE_ID)-5)

unique(df$TIMEPOINT)
df[df$TIMEPOINT == 'ND'] <- 'Newly Diagnosed'

df$TIMEPOINT[df$TIMEPOINT == 'ND'] <- 'Newly Diagnosed'

unique(df$TISSUE)
df$TISSUE[df$TISSUE == 'BM'] <- 'Bone Marrow'
df$TISSUE[df$TISSUE == 'PB'] <- 'Peripheral Blood'

df$CANCER_TYPE_DETAILED <- ifelse(df$CANCER_TYPE == 'MM', 'Multiple Myeloma', df$CANCER_TYPE_DETAILED)
df$CANCER_TYPE_DETAILED <- ifelse(df$CANCER_TYPE == 'PCL', 'Plasma Cell Leukemia', df$CANCER_TYPE_DETAILED)

df$CANCER_TYPE_DETAILED <- ifelse(df$CANCER_TYPE == 'MM', 'Multiple Myeloma', )

df$CANCER_TYPE_DETAILED <- ifelse(df$CANCER_TYPE == 'MM', 'Multiple Myeloma',
                                  ifelse(df$CANCER_TYPE == 'PCL', 'Plasma Cell Leukemia',
                                         df$CANCER_TYPE))

df <- dplyr::select(df, c('PATIENT_ID', 'SAMPLE_ID', 'TISSUE', 'CANCER_TYPE',
                          'CANCER_TYPE_DETAILED', 'TIMEPOINT', 'TRANSLOCATION'))
df_old <- read.csv('data_clinical_sample.txt', sep="\t", header=TRUE)

df_old <- read.csv('data_clinical_sample.txt', sep="\t", header=FALSE, comment.char = "~")

library(purrr)
library(dplyr)
library(stringr)

setwd("C:/Users/abhmalat/OneDrive - Indiana University/cbio_MMRF")

df <- read.csv('raw_files/MMRF_IA15_summary.csv', sep=",", header=FALSE)

names(df) <- c("SAMPLE_ID","Patient","TISSUE","TIMEPOINT","CANCER_TYPE","t(4;14) cons",
               "t(11;14) cons","t(6;14) cons","t(14;16) cons","t(14;20) cons","t(8;14) cons",
               "None","Hyperdiploid","MYC","Mutational Signature")

head(df,10)

df <- df[2:nrow(df), ]

df[df == '#N/A'] <- 'NA'
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

df_old <- read.csv('data_clinical_sample.txt', sep="\t", header=FALSE, comment.char = "~")

df$CANCER_TYPE_DETAILED <- ifelse(df$CANCER_TYPE == 'MM', 'Multiple Myeloma',
                                  ifelse(df$CANCER_TYPE == 'PCL', 'Plasma Cell Leukemia',
                                         df$CANCER_TYPE))

df <- dplyr::select(df, c('PATIENT_ID', 'SAMPLE_ID', 'TISSUE', 'CANCER_TYPE',
                          'CANCER_TYPE_DETAILED', 'TIMEPOINT', 'TRANSLOCATION'))

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

df[df == '#N/A'] <- 'NA'
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

df_old <- read.csv('data_clinical_sample.txt', sep="\t", header=FALSE, comment.char = "%")


df_final <- rbind(dplyr::select(df, c('PATIENT_ID', 'SAMPLE_ID', 'TISSUE', 'CANCER_TYPE',
                                      'CANCER_TYPE_DETAILED', 'TIMEPOINT', 'TRANSLOCATION',
                                      df_old)
df_final <- rbind(dplyr::select(df, c('PATIENT_ID', 'SAMPLE_ID', 'TISSUE', 'CANCER_TYPE',
                                      'CANCER_TYPE_DETAILED', 'TIMEPOINT', 'TRANSLOCATION',
                                      df_old))
df_final <- rbind(dplyr::select(df, c('PATIENT_ID', 'SAMPLE_ID', 'TISSUE', 'CANCER_TYPE',
                                      'CANCER_TYPE_DETAILED', 'TIMEPOINT', 'TRANSLOCATION',
                                      df_old)))
colnames(df_old) <- df_old(5, )

colnames(df_old) <- df_old[5, ]

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

df[df == '#N/A'] <- 'NA'
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

df_final <- rbind(dplyr::select(df, c('PATIENT_ID', 'SAMPLE_ID', 'TISSUE', 'CANCER_TYPE',
                                      'CANCER_TYPE_DETAILED', 'TIMEPOINT', 'TRANSLOCATION',
                                      df_old)))

df_final <- rbind(dplyr::select(df, c('PATIENT_ID', 'SAMPLE_ID', 'TISSUE', 'CANCER_TYPE',
                                      'CANCER_TYPE_DETAILED', 'TIMEPOINT', 'TRANSLOCATION'),
                                      df_old))
df %>% dplyr::select(c('PATIENT_ID', 'SAMPLE_ID', 'TISSUE', 'CANCER_TYPE',
                                      'CANCER_TYPE_DETAILED', 'TIMEPOINT', 'TRANSLOCATION')
df_final <- rbind(df %>% select(c('PATIENT_ID', 'SAMPLE_ID', 'TISSUE', 'CANCER_TYPE',
                                      'CANCER_TYPE_DETAILED', 'TIMEPOINT', 'TRANSLOCATION')),
                                      df_old) %>% order('SAMPLE_ID')
df_final <- rbind(df %>% select(c('PATIENT_ID', 'SAMPLE_ID', 'TISSUE', 'CANCER_TYPE',
                                      'CANCER_TYPE_DETAILED', 'TIMEPOINT', 'TRANSLOCATION')),
                                      df_old)
df_final <- df_old[1:5, ]

df_final <- rbind(df %>% select(c('PATIENT_ID', 'SAMPLE_ID', 'TISSUE', 'CANCER_TYPE',
                                      'CANCER_TYPE_DETAILED', 'TIMEPOINT', 'TRANSLOCATION')),
                                      df_old[6:nrow(df_old),])

df_final <- df_old[1:5, ]
df_final <- rbind(df_final, df %>% select(c('PATIENT_ID', 'SAMPLE_ID', 'TISSUE', 'CANCER_TYPE',
                                      'CANCER_TYPE_DETAILED', 'TIMEPOINT', 'TRANSLOCATION')),
                                      df_old[6:nrow(df_old),])
df_final <- unique(df_final, by = c('SAMPLE_ID', 'PATIENT_ID'))

df_final[df_final == ''] <- NA


df_final <- left_join(df_final, df[, c('SAMPLE_ID', 'PATIENT_ID', 'HYPERDIPLOID','MYC','MUTATIONAL_SIGNATURE')],
                      by = c('SAMPLE_ID', 'PATIENT_ID'))
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
outputFile <- "data_clinical_sample_1203.txt"

df_final[1, c('HYPERDIPLOID','MYC','MUTATIONAL_SIGNATURE')]

df_final[1, c('HYPERDIPLOID','MYC','MUTATIONAL_SIGNATURE')] <- c('Hyperdiploid','MYC','Mutational Signature')

df_final[2, c('HYPERDIPLOID','MYC','MUTATIONAL_SIGNATURE')] <- c('Hyperdiploid','MYC','Mutational Signature')
df_final[3, c('HYPERDIPLOID','MYC','MUTATIONAL_SIGNATURE')] <- c('STRING', 'STRING', 'STRING')
df_final[4, c('HYPERDIPLOID','MYC','MUTATIONAL_SIGNATURE')] <- c(8,9,10)
df_final[5, c('HYPERDIPLOID','MYC','MUTATIONAL_SIGNATURE')] <- c('HYPERDIPLOID','MYC','MUTATIONAL_SIGNATURE')

write.table(df_final, file=outputFile, sep = "\t", col.names = FALSE,
            row.names=FALSE, quote = FALSE, append = FALSE)
print('Samples completed')

metaFile <- file("meta_patient_sample.txt")
writeLines(c(
  "cancer_study_identifier: mmrf_2020",
  "genetic_alteration_type: CLINICAL",
  "datatype: SAMPLE_ATTRIBUTES",
  paste("data_filename: ", outputFile)
), metaFile
)
close(metaFile)

metaFile <- file("meta_clinical_sample.txt")
writeLines(c(
  "cancer_study_identifier: mmrf_2020",
  "genetic_alteration_type: CLINICAL",
  "datatype: SAMPLE_ATTRIBUTES",
  paste("data_filename: ", outputFile)
), metaFile
)
close(metaFile)
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


df <- read.csv('data_timeline_lab_test_ALL_MM.txt', sep="\t", header=TRUE, quote = "")
write.table(df, file='data_timeline_lab_test_ALL_MM_formatted.txt', sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")


df_treatment <- read.csv('data_timeline_treatment_ALL_MM.txt', sep="\t", header=TRUE)
write.table(df_treatment, file='data_timeline_treatment_ALL_MM_formatted.txt', sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")


library(purrr)
library(org.Hs.eg.db)
library(annotate)
library(dplyr)
library(tidyr)
library(biomaRt)

setwd('c:/Users/abhmalat/OneDrive - Indiana University/cbio_MMRF')

df <- read.csv('somatic_mutation_files/MMRF_CoMMpass_IA16a_All_Canonical_NS_Variants.txt', sep="\t", header=TRUE)


# SEPARATE ANN....EFFECT field into individual columns
df <- separate(data=df, col=ANN....EFFECT, into = c("ANN_EFFECT1","ANN_EFFECT2", "ANN_EFFECT3", "ANN_EFFECT4"), sep="&", remove=FALSE, fill= "right")

levels(df$ANN....BIOTYPE)

ensembl <- useMart(host="apr2020.archive.ensembl.org",
                   biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')


Entrez_Gene_Ids <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), filters = 'ensembl_gene_id', values=df$ANN....GENEID, ensembl)

Entrez_Gene_Ids[,"ensembl_gene_id"] <- as.factor(Entrez_Gene_Ids[,"ensembl_gene_id"])

Entrez_Gene_Ids[,"entrezgene_id"] <- as.factor(Entrez_Gene_Ids[,"entrezgene_id"])

foo <- unique(arrange(Entrez_Gene_Ids, ensembl_gene_id, entrezgene_id))

dfEntrez <- left_join(df, foo, by=c("ANN....GENEID" = "ensembl_gene_id") )


dfVariants <- data.frame("ANN_EFFECT1" = c("disruptive_inframe_insertion", "inframe_insertion", "disruptive_inframe_deletion", "inframe_deletion", "missense_variant", "stop_gained", "synonymous_variant", "stop_retained_variant", "splice_donor_variant", "splice_acceptor_variant", "start_lost", "initiator_codon_variant", "stop_lost", "3_prime_UTR_variant", "downstream_gene_variant", "5_prime_UTR_variant", "5_prime_UTR_premature_start_codon_gain_variant", "5_prime_UTR_truncation", "upstream_gene_variant", "non_coding_exon_variant", "exon_loss_variant", "intron_variant", "intragenic_variant", "splice_region_variant", "frameshift_variant"),
                         "Variant_Classification" = c("In_Frame_Ins", "In_Frame_Ins", "In_Frame_Del", "In_Frame_Del", "Missense_Mutation", "Nonsense_Mutation", "Silent", "Silent", "Splice_Site", "Splice_Site", "Translation_Start_Site", "Translation_Start_Site", "Nonstop_Mutation", "3'UTR", "3'Flank", "5'UTR", "5'UTR", "5'UTR", "5'Flank", "Intron", "Splice_Site", "Intron", "Intron", "Splice_Region", "Unknown"))
head(dfVariants, 50)

dfVarEntrez <- left_join(dfEntrez, dfVariants, by="ANN_EFFECT1")

dfVarEntrez2 <- dfVarEntrez %>% dplyr::select( ANN....GENE, entrezgene_id, X.CHROM, POS, Variant_Classification, REF, ALT, ID, Sample, ANN....HGVS_P, ANN....AA_POS )


names(dfVarEntrez2) <- c("Hugo_Symbol", "Entrez_Gene_Id", "Chromosome", "Start_Position", "Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele1", "dbSNP_RS", "Tumor_Sample_Barcode", "HGVSp_Short", "Protein_Pos")


aaa0<-as.matrix(dfVarEntrez2)

ddd<-nchar(aaa0[,"Reference_Allele"])-nchar(aaa0[,"Tumor_Seq_Allele1"])
ddd[which(ddd<0)]<-0
ddd1<-nchar(aaa0[,"Reference_Allele"])-nchar(aaa0[,"Tumor_Seq_Allele1"])
Start_Position<-as.numeric(aaa0[,"Start_Position"])
End_Position<-as.numeric(aaa0[,"Start_Position"])+ddd

eee<-rep("",length(ddd1))
eee[which(ddd1==0)]<-"SNP"
eee[which(ddd1>0)]<-"DEL"
eee[which(ddd1<0)]<-"INS"


dfVarEntrez2$End_Position = End_Position

dfVarEntrez2$Variant_Type = eee

dfVarEntrez2$NCBI_Build='GRCh37'

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

library(purrr)
library(dplyr)
library(org.Hs.eg.db)
library(annotate)

setwd('C:/Users/abhmalat/OneDrive - Indiana University/cbio_MM/')

# Fix quotes for height
df <- read.csv('data_timeline_lab_test_ALL_MM.txt', sep="\t", header=TRUE, quote = "")
write.table(df, file='data_timeline_lab_test_ALL_MM_formatted.txt', sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")


df_treatment <- read.csv('data_timeline_treatment_ALL_MM.txt', sep="\t", header=TRUE)
write.table(df_treatment, file='data_timeline_treatment_ALL_MM_formatted.txt', sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")



setwd('C:/Users/abhmalat/OneDrive - Indiana University/cbio_MM/')

# Fix quotes for height
df_labtest <- read.csv('data_timeline_lab_test_ALL_MM.txt', sep="\t", header=TRUE, quote = "")
df_labtest_extra_header <- read.csv('data_extra_header_timeline_lab_test_ALL_MM.txt', sep="\t", header=TRUE, quote = "")

lt_headers <- cbind(df_labtest[1,1:6], lt_headers[1,])
lt_headers <- cbind(df_labtest[1,1:6], df_labtest_extra_header[1,])


# Fix quotes for height
labtest_raw <- read.csv('data_timeline_lab_test_ALL_MM.txt', sep="\t", header=TRUE, quote = "")
lt_header <- read.csv('data_extra_header_timeline_lab_test_ALL_MM.txt', sep="\t", header=TRUE, quote = "")

lt_colnames <- cbind(labtest_raw[1, 1:6], lt_header[1,])

labtest_raw <- read.csv('data_timeline_lab_test_ALL_MM.txt', sep="\t", header=FALSE, quote = "")
lt_header <- read.csv('data_extra_header_timeline_lab_test_ALL_MM.txt', sep="\t", header=FALSE, quote = "")

lt_colnames <- cbind(labtest_raw[1, 1:6], lt_header[1,])
labtests <- cbind(lt_colnames, labtest_raw)

labtest <- read.csv("data_timeline_lab_test_ALL_PST.txt", sep = "\t", header = FALSE,
                    col.names = toupper(lt_colnames), na.strings = c("N/A", "", "unavailable"))
# Fix quotes for height
labtest_raw <- read.csv('data_timeline_lab_test_ALL_MM.txt', sep="\t", header=FALSE, quote = "")
lt_header <- read.csv('data_extra_header_timeline_lab_test_ALL_MM.txt', sep="\t", header=FALSE, quote = "")

lt_colnames <- cbind(labtest_raw[1, 1:6], lt_header[1,])
labtest[!is.na(labtest$grade),]
labtest[,!is.na(labtest$grade)]
write.table(labtest, file='data_timeline_lab_test_ALL_MM_formatted.txt', sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")

library(purrr)
library(dplyr)
library(org.Hs.eg.db)
library(annotate)

setwd('C:/Users/abhmalat/OneDrive - Indiana University/cbio_MM/')

# Fix quotes for height
labtest_raw <- read.csv('data_timeline_lab_test_ALL_MM.txt', sep="\t", header=FALSE, quote = "")
lt_header <- read.csv('data_extra_header_timeline_lab_test_ALL_MM.txt', sep="\t", header=FALSE, quote = "")

lt_colnames <- cbind(labtest_raw[1, 1:7], lt_header[1,])
labtest <- read.csv("data_timeline_lab_test_ALL_MM.txt", sep = "\t", header = FALSE,
                    col.names = toupper(lt_colnames), na.strings = c("N/A", "", "unavailable"))
write.table(labtest, file='data_timeline_lab_test_formatted.txt', sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")


df_treatment <- read.csv('data_timeline_treatment_ALL_MM.txt', sep="\t", header=TRUE)
write.table(df_treatment, file='data_timeline_treatment_formatted.txt', sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")


labtest[!is.na(labtest$GRADE), ]
ltnew <- labtest[2:nrow(labtest), ]

ltnew[35624,]
ltnew[35625,]
ltnew[35623,]
ltnew[35623,"RESULT"]
labtest_raw[35623,"RESULT"]
labtest_raw[35622:35627,"RESULT"]
labtest_raw <- read.csv('data_timeline_lab_test_ALL_MM.txt', sep="\t", header=FALSE, quote = "")

lt_colnames <- cbind(labtest_raw[1, 1:7], lt_header[1,])

labtest <- read.csv("data_timeline_lab_test_ALL_MM.txt", sep = "\t", header = FALSE,
                    col.names = toupper(lt_colnames), na.strings = c("N/A", "", "unavailable"))
labtest <- read.csv("data_timeline_lab_test_ALL_MM.txt", sep = "\t", header = FALSE,
                    col.names = toupper(lt_colnames), na.strings = c("N/A", "", "unavailable"), quote = TRUE)
labtest <- read.csv("data_timeline_lab_test_ALL_MM.txt", sep = "\t", header = FALSE,
                    col.names = toupper(lt_colnames), na.strings = c("N/A", "", "unavailable"), quote = '"')
labtest <- read.csv("data_timeline_lab_test_ALL_MM.txt", sep = "\t", header = FALSE,
                    col.names = toupper(lt_colnames), na.strings = c("N/A", "", "unavailable"), quote = '~')
labtest[35622:35627,"RESULT"]
write.table(labtest[2:nrow(labtest), ], file='data_timeline_lab_test_formatted.txt', sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")

df_treatment1 <- unique(df_treatment)
write.table(unique(df_treatment), file='data_timeline_treatment_formatted.txt', sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")

labtest1 <- unique(labtest[2:nrow(labtest), ])
(write.table(unique(labtest[2:nrow(labtest), ]), file='data_timeline_lab_test_formatted.txt', sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")
write.table(unique(labtest[2:nrow(labtest), ]), file='data_timeline_lab_test_formatted.txt', sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")
df_treatment1[(df_treatment$PATIENT_ID == '0661-515'), ]
df_treatment1[(df_treatment$PATIENT_ID == '0661-515' & df_treatment1$AGENT == 'GABAPENTIN CAP'), ]
df_treatment1[df_treatment$PATIENT_ID == '0661-515' & df_treatment1$AGENT == 'GABAPENTIN CAP', ]
df_treatment1[df_treatment1$AGENT == 'GABAPENTIN CAP', ]
library(stringr)
dosage <- str_split(df_treatment1$DOSAGE, ',')
dosage <- as.data.frame(str_split(df_treatment1$DOSAGE, ','))
library(data.table)
dosage <- as.data.table(str_split(df_treatment1$DOSAGE, ','))
dosage <- df_treatment1$DOSAGE
dosage <- as.data.frame.matrix(df_treatment1$DOSAGE)
dosage <- as.data.frame(df_treatment1$DOSAGE)
colnames(dosage) <- "V1"
library(dplyr)
library(tidyr)
labtest <- unique(labtest)
write.table(unique(labtest[2:nrow(labtest), c("PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "TEST", "RESULT", "EVENT_DATE", "RISK_EVALUATION", "AFFECTED_CELLS", "TOTAL_CELLS", "GRADE")]), file='data_timeline_lab_test_formatted.txt', sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")

a0661 <- unique(labtest[2:nrow(labtest), c("PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "TEST", "RESULT", "EVENT_DATE", "RISK_EVALUATION", "AFFECTED_CELLS", "TOTAL_CELLS", "GRADE")]) %>% dplyr::filter("PATIENT_ID" = "0661-515" & "TEST" = "Platelets [#/volume] in Blood k/cumm")
a0661 <- unique(labtest[2:nrow(labtest), c("PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "TEST", "RESULT", "EVENT_DATE", "RISK_EVALUATION", "AFFECTED_CELLS", "TOTAL_CELLS", "GRADE")]) %>% dplyr::filter(c("PATIENT_ID" = "0661-515", "TEST" = "Platelets [#/volume] in Blood k/cumm"))
a0661 <- unique(labtest[2:nrow(labtest), c("PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "TEST", "RESULT")]) %>% dplyr::filter(c("PATIENT_ID" == "0661-515", "TEST" == "Platelets [#/volume] in Blood k/cumm"))
a0661 <- labtest %>% dplyr::filter(c("PATIENT_ID" == "0661-515", "TEST" == "Platelets [#/volume] in Blood k/cumm"))
a0661 <- labtest %>% dplyr::filter("PATIENT_ID" == "0661-515" & "TEST" == "Platelets [#/volume] in Blood k/cumm")
labtest %>% dplyr::filter("PATIENT_ID" == "0661-515" & "TEST" == "Platelets [#/volume] in Blood k/cumm")
labtest %>% dplyr::filter("PATIENT_ID" == "0661-515")
labtest %>% dplyr::filter("PATIENT_ID" = "0661-515")
labtest %>% dplyr::filter(PATIENT_ID == "0661-515")
a0661 <- labtest %>% dplyr::filter(PATIENT_ID == "0661-515" & TEST == "Platelets [#/volume] in Blood k/cumm")
write.table(unique(labtest[2:nrow(labtest), c("PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "TEST", "RESULT", "EVENT_DATE", "RISK_EVALUATION", "AFFECTED_CELLS", "TOTAL_CELLS", "GRADE")]) %>% dplyr::filter(PATIENT_ID == "0661-515" & TEST == "Platelets [#/volume] in Blood k/cumm"), file='data_timeline_lab_test_formattedtest.txt', sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")

library(purrr)
library(dplyr)


setwd('C:/Users/abhmalat/OneDrive - Indiana University/cbio_MM/')
df_treatment <- read.csv('data_timeline_treatment_ALL_MM.txt', sep="\t", header=TRUE)

trt_u <- unique(df_treatment)

dup <- subset(df_treatment, !(df_treatment %in% trt_u))
dup <- subset(trt_u, !(df_treatment %in% trt_u))
dup <- !(df_treatment %in% trt_u)
dup <- setdiff(df_treatment, trt_u)
dup <- setdiff(trt_u, df_treatment)
dup <- anti_join(df_treatment, trt_u)
dup <- duplicated(df_treatment)
dup <- df_treatment[dup]
dup <- df_treatment[dup, ]
write.table(dup, "duplicated_treatments.txt", sep="\t", col.names = TRUE, row.names = FALSE, quote = "")
write.table(dup, "duplicated_treatments.txt", sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
setwd("c:/Users/abhmalat/OneDrive - Indiana University/cbio_MMRF")
list.files()
sampleinfo <- read.csv("0661 sample info.xlsx", sep = ",", header = TRUE)
sampleinfo <- read.csv("extra_sample_info.csv", sep = ",", header = TRUE)
samples <- read.csv("data_clinical_sample_ALL_MM.txt", sep = "\t", header = FALSE)
sampleinfo <- read.csv("extra_sample_info.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
sampleinfo <- read.csv("extra_sample_info.csv", sep = ",", header = TRUE, stringsAsFactors = TRUE)
sampleinfo$Sample.name <- as.character(sampleinfo$Sample.name)
sampleinfo <- read.table("extra_sample_info.csv", sep = ",", header = TRUE, stringsAsFactors = TRUE)
sampleinfo <- read.table("extra_sample_info.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
sampleinfo <- read.table("extra_sample_info.csv", sep = ",", header = TRUE, colClasses = "character")
setdiff(sampleinfo$Sample.name, samples$V2)
samples$V2
library(dplyr)
library(stringr)
library(data.table)
library(maftools)

setwd("c:/Users/abhmalat/OneDrive - Indiana University/cbio_MM")
smpl <- read.csv("data_clinical_sample_ALL_MM.txt", sep="\t", header=FALSE)
smpl <- read.csv("data_clinical_sample_ALL_MM.txt", sep="\t", header=FALSE, na="")
intermediate <- smpl[6:nrow(smpl),]
colnames(intermediate) <- smpl[5,]
som <- read.csv("data/somatic/vcf_metadata_mm_somatic1612815273.csv", sep=",", header = TRUE)
setdiff(som$tumor_id, intermediate$SAMPLE_ID)
setdiff(som$normal_id, intermediate$SAMPLE_ID)
unique(setdiff(som$normal_id, intermediate$SAMPLE_ID))
germ$sample <- lapply(germ$normal_id, function(x) {x[1:str_locate(x, "-")]})
str_locate(germ$normal_id, "-")
str_locate(germ$normal_id, "-")[1]
germ$sample <- lapply(germ$normal_id, function(x) {x[1:10]})
lapply(germ$normal_id, function(x) {x[1:10]})
lapply(germ$normal_id, function(x) {x[1:10]}[1])
germ$sample <- lapply(germ$normal_id, function(x) {x[1:10]}[1])
str_extract(germ$normal_id, "-")
str_extract(germ$normal_id, "*-*")
str_extract(germ$normal_id, "*.-.*")
str_extract(germ$normal_id, "*.-")
germ <- read.table("data/germline/vcf_metadata_germline.csv", sep=",", header = TRUE, colClasses = "character")
?lapply
lapply(germ$normal_id, function (x) {substr(x,1,str_locate(x,"-"))})
lapply(germ$normal_id, function (x) {substr(x,1,str_locate(x,"-")-1)})
germ$sample <- lapply(germ$normal_id, function (x) {substr(x,1,str_locate(x,"-")-1)})
walker <- read.csv("0661 sample info.xlsx", sep=",", header = TRUE)
walker <- read.csv("0661 sample info.xlsx", sep=",")
walker <- read.csv("extra_sample_info.csv", sep=",", header = TRUE)
walker <- read.table("extra_sample_info.csv", sep=",", header = TRUE, colClasses = "character")
setdiff(intermediate$SAMPLE_ID, walker$Sample.name)
setdiff(walker$Sample.name, intermediate$SAMPLE_ID)
setdiff(walker$Sample.name, som$tumor_id)
setdiff(walker$Sample.name, som$normal_id)
rm(walker)
setdiff(smpl$V2, germ$sample_id)
setdiff(intermediate$SAMPLE_ID, germ$sample_id)
setdiff(germ$sample_id, intermediate$SAMPLE_ID)
colnames(germ)[2]
colnames(germ)[2] <- "VENDOR"
germ$TUMOR_TISSUE_SITE <- NA
germ$TUMOR_TYPE <- NA
germ$SAMPLE_TYPE <- NA
colnames(germ)
germ_fixed <- germ %>% dplyr::select(c("patient_id", "sample_id", "SAMPLE_TYPE", "TUMOR_TISSUE_SITE", "TUMOR_TYPE", "VENDOR"))
colnames(germ_fixed)[1] <- "PATIENT_ID"
colnames(germ_fixed)[2] <- "SAMPLE_ID"
prefinal <- rbind(intermediate, germ_fixed)
prefinal <- unique(rbind(intermediate, germ_fixed))
som$PATIENT_ID <- lapply(som$filename, function (x) {substr(x,1,str_locate(x,"-")-1)})
som$PATIENT_ID <- lapply(som$filename, function (x) {substr(x,1,str_locate(x,"_")-1)})
som_fixed <- som %>% dplyr::select(c(PATIENT_ID, tumor_id))
som_fixed$TUMOR_TISSUE_SITE <- NA
som_fixed$SAMPLE_TYPE <- NA
som_fixed$TUMOR_TYPE <- NA
som_fixed$VENDOR <- "CMG"
colnames(som_fixed)
colnames(som_fixed)[2] <- "SAMPLE_ID"
intermediate <- unique(rbind(intermediate, som_fixed))
duplicated(intermediate[,1:2])
intermediate[duplicated(intermediate[,1:2]),]
setdiff(intermediate, intermediate[duplicated(intermediate[,1:2]),])
final <- setdiff(intermediate, intermediate[duplicated(intermediate[,1:2]),])
colnames(smpl) <- smpl[5,]
realfinal <- rbind(smpl[1:5], final)
smpl[1:5]
smpl[1:5,]
realfinal <- rbind(smpl[1:5,], final)
write.table(realfinal, "data_clinical_sample_formatted.txt", col.names = FALSE, row.names = FALSE, sep="\t", quote = FALSE, append = FALSE, na = "")
write.table(as.data.frame(realfinal), "data_clinical_sample_formatted.txt", col.names = FALSE, row.names = FALSE, sep="\t", quote = FALSE, append = FALSE, na = "")
fnl <- apply(realfinal, 2, as.character)
write.csv(fnl, "data_clinical_sample_formatted.txt", col.names = FALSE, row.names = FALSE, sep="\t", quote = FALSE, append = FALSE, na = "")
write.table(fnl, "data_clinical_sample_formatted.txt", col.names = FALSE, row.names = FALSE, sep="\t", quote = FALSE, append = FALSE, na = "")
sv <- som %>% dplyr::filter(vcf_type=='somaticSV')
snvs <- som %>% dplyr::filter(vcf_type=='snvs')
indels <- som %>% dplyr::filter(vcf_type=='indels')
mutationCL <- ("case_lists/cases_sv.txt")
getwd()
writeLines(c(
  "cancer_study_identifier: PHI_MM_2020",
  "stable_id: PHI_MM_2020_sv",
  "case_list_name: SomaticSV",
  "case_list_description: Somatic Structural Variants",
  paste("case_list_ids: ", paste(unique(sv$tumor_id), collapse = '\t'))
  ), f
)
mutationCL <- ("case_lists/cases_snvs.txt")
writeLines(c(
  "cancer_study_identifier: PHI_MM_2020",
  "stable_id: PHI_MM_2020_snvs",
  "case_list_name: SNVS",
  "case_list_description: Single Nucleotide Variants",
  paste("case_list_ids: ", paste(unique(snvs$tumor_id), collapse = '\t'))
  ), f
)
mutationCL <- ("case_lists/cases_indels.txt")
writeLines(c(
  "cancer_study_identifier: PHI_MM_2020",
  "stable_id: PHI_MM_2020_indels",
  "case_list_name: INDELS",
  "case_list_description: Insertion/Deletions",
  paste("case_list_ids: ", paste(unique(indels$tumor_id), collapse = '\t'))
  ), f
)
mutationCL <- ("case_lists/cases_germline.txt")
f <- file(mutationCL)
writeLines(c(
  "cancer_study_identifier: PHI_MM_2020",
  "stable_id: PHI_MM_2020_germline",
  "case_list_name: Germline",
  "case_list_description: Germline Samples",
  paste("case_list_ids: ", paste(unique(germ$sample_id), collapse = '\t'))
  ), f
)