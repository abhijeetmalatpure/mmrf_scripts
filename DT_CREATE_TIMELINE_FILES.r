library(purrr)
library(dplyr)
library(org.Hs.eg.db)
library(annotate)

setwd('f:/mmrf_data')

# TREATMENT
df <- read.csv('clinical_flat_files/MMRF_CoMMpass_IA11_STAND_ALONE_TREATMENT_REGIMEN.csv', sep=",", header=TRUE)

head(df)
summary(df$stopday)
class(df$stopday)

df$EVENT_TYPE = "TREATMENT"
df$TREATMENT_TYPE = "Medical Therapy"
df$SUBTYPE = "Chemotherapy"

df["startday"][is.na(df["startday"])] <- -1
df["stopday"][is.na(df["stopday"])] <- -1

dfFinal <- df %>% dplyr::select(public_id, startday, stopday, EVENT_TYPE, TREATMENT_TYPE, SUBTYPE, MMTX_THERAPY)

names(dfFinal) <- c("PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "TREATMENT_TYPE", "SUBTYPE", "AGENT")

head(dfFinal,100)

outputFile <- "DT_R/cbio_files/data_timeline_treatment.txt"

file.remove(outputFile)
file.create(outputFile)

write.table(dfFinal, file=outputFile, sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="") 


# STATUS
df <- read.csv('clinical_flat_files/MMRF_CoMMpass_IA11_STAND_ALONE_TRTRESP.csv', sep=",", header=TRUE)

head(df)
summary(df)

df$EVENT_TYPE = "STATUS"

df["therstdy"][is.na(df["therstdy"])] <- -1
df["therendy"][is.na(df["therendy"])] <- -1

dfFinal <- df %>% dplyr::select(public_id, therstdy, therendy, EVENT_TYPE, bestresp)

names(dfFinal) <- c("PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "STATUS")

head(dfFinal,100)

outputFile <- "DT_R/cbio_files/data_timeline_status.txt"

file.remove(outputFile)
file.create(outputFile)

write.table(dfFinal, file=outputFile, sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, na="") 



