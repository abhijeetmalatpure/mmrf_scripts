library(purrr)
library(dplyr)


setwd('C:/Users/abhmalat/OneDrive - Indiana University/cbio_MM/')

# Fix quotes for height
labtest_raw <- read.csv('data_timeline_lab_test_ALL_MM.txt', sep="\t", header=FALSE, quote = "")
lt_header <- read.csv('data_extra_header_timeline_lab_test_ALL_MM.txt', sep="\t", header=FALSE, quote = "")

lt_colnames <- cbind(labtest_raw[1, 1:7], lt_header[1,])
labtest <- read.csv("data_timeline_lab_test_ALL_MM.txt", sep = "\t", header = FALSE,
                    col.names = toupper(lt_colnames), na.strings = c("N/A", "", "unavailable"), quote = '~')
ltnew <- labtest[2:nrow(labtest), ]
write.table(unique(labtest[2:nrow(labtest), ]), file='data_timeline_lab_test_formattedtest.txt',
            sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")


df_treatment <- read.csv('data_timeline_treatment_ALL_MM.txt', sep="\t", header=TRUE)
trt_u <- unique(df_treatment)
write.table(unique(df_treatment), file='data_timeline_treatment_formatted.txt',
            sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="")

