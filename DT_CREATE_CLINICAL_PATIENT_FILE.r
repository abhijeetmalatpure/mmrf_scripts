library(dplyr)
library(tidyr)

setwd("c:/Users/abhmalat/OneDrive - Indiana University/RI cBioPortal/PEDS_brrhelm")

per_patient <- read.csv('', header=TRUE)

head(per_patient)

colnames(per_patient)

pat_cols <- c('PUBLIC_ID', 'D_PT_DIDPATIENTCOM', 'D_PT_PRIMARYREASON', 'D_PT_CAUSEOFDEATH', 'D_PT_deathdy', 'D_PT_trtstdy', 'D_PT_sdeathdy', "D_PT_mmstatus", "D_PT_mmstatus1", "D_PT_mmstatus2", "D_PT_mmstatus3", 'D_PT_respdur', 'DEMOG_GENDER', 'D_PT_age', 'D_PT_race', 'D_PT_ethnic', 'D_PT_gender')

per_pat_keep <- per_patient[, pat_cols]

head(per_pat_keep)

colnames(per_pat_keep)

summary(per_pat_keep)

nrow(per_pat_keep)

summary(per_pat_keep)

mmrf_pat_keep_names <- c('PUBLIC_ID', 'D_PT_DIDPATIENTCOM', 'D_PT_PRIMARYREASON', 'D_PT_deathdy', 'D_PT_sdeathdy', 'D_PT_trtstdy', 'D_PT_respdur', 'PUBLIC_ID', 'DEMOG_GENDER', 'D_PT_age', 'D_PT_race', 'D_PT_ethnic', 'D_PT_gender')

cbio_pat_keep_names <- c('PATIENT_ID', 'OS_STATUS', 'OS_STATUS_SUP', 'OS_MONTHS_ALT', 'OS_MONTHS', 'DFS_STATUS', 'DFS_MONTHS', 'PATIENT_DISPLAY_NAME', 'GENDER_PART', 'AGE', 'RACE_RAW', 'ETHNICITY_RAW', 'GENDER_RAW')

racedf <- data.frame(RACE_RAW = c(1, 2, 3, 4, 5, 6), RACE = c("white", "black/african", "american indian", "asian", "native hawaiian", "other"))

genderdf <- data.frame(GENDER_RAW = c(1, 2), GENDER=c("Male", "Female"))

ethnicitydf <- data.frame(ETHNICITY_RAW = c(1, 2, 3), ETHNICITY=c("Hispanic", "Non-hispanic", "Other"))

cbio_pat <- per_pat_keep[, mmrf_pat_keep_names]

head(cbio_pat)

names(cbio_pat) <- cbio_pat_keep_names

head(cbio_pat)

cbio_pat <- left_join(cbio_pat, racedf, by="RACE_RAW") %>% dplyr::select(-RACE_RAW)

cbio_pat <- left_join(cbio_pat, genderdf, by="GENDER_RAW") %>% dplyr::select(-GENDER_RAW, -GENDER_PART)

cbio_pat <- left_join(cbio_pat, ethnicitydf, by="ETHNICITY_RAW") %>% dplyr::select(-ETHNICITY_RAW)


summary(cbio_pat[,c('GENDER', 'RACE', 'ETHNICITY')])


head(cbio_pat)

cbio_pat$OS_STATUS <- with(cbio_pat, ifelse(OS_STATUS == 'No' & OS_STATUS_SUP == 'Death', 'DECEASED', 'LIVING'))
cbio_pat$DFS_STATUS <- with(cbio_pat, ifelse(DFS_STATUS == 1, 'DiseaseFree', 'Recurred/Progressed'))

cbio_pat$OS_MONTHS <- with(cbio_pat,
                         ifelse(
                             is.na(OS_MONTHS) & OS_STATUS == 'LIVING',
                                 trunc(max(per_pat_keep$D_PT_deathdy,na.rm = TRUE)/30),
                                 ifelse(
                                     is.na(OS_MONTHS),
                                         trunc(OS_MONTHS_ALT/30),
                                         trunc(OS_MONTHS/30)
                                     )
                             )
                         )



cbio_pat$DFS_MONTHS <- with(cbio_pat, ifelse(is.na(DFS_MONTHS), 0, trunc(DFS_MONTHS/30)))

output_cols <- c('PATIENT_ID', 'OS_STATUS', 'OS_MONTHS', 'DFS_STATUS', 'DFS_MONTHS', 'PATIENT_DISPLAY_NAME', 'GENDER', 'AGE', 'RACE', 'ETHNICITY')

head(cbio_pat[,output_cols])

summary(cbio_pat[,output_cols])

head(cbio_pat[,output_cols], 100)

outputFile <- "DT_R/cbio_files/data_clinical_patient.txt"

file.remove(outputFile)
file.create(outputFile)
f <- file(outputFile)
writeLines(c("#Patient Identifier\tOverall Survival Status\tOverall Survival (Months)\tDisease Free Status\tDisease Free (Months)\tPatient Name\tGender\tAge\tRace\tEthnicity",
             "#ID of Patient in Study\tOverall Patient Survival Status\tOverall Patient Survival Months\tDisease Free Status from initial treatment\tDisease Free Months since initial treatment\tName or ID of Patient\tGender of Patient\tAge of Patient when first diagnosed\tRace\tEthnicity",
             "#STRING\tSTRING\tNUMBER\tSTRING\tNUMBER\tSTRING\tSTRING\tNUMBER\tSTRING\tSTRING",
             "#1\t2\t3\t4\t5\t6\t7\t8\t9\t10",
             "PATIENT_ID\tOS_STATUS\tOS_MONTHS\tDFS_STATUS\tDFS_MONTHS\tPATIENT_DISPLAY_NAME\tGENDER\tAGE\tRACE\tETHNICITY"), f)
close(f)

write.table(cbio_pat[,output_cols], file=outputFile, sep = "\t", col.names = FALSE, row.names=FALSE, quote = FALSE, append = TRUE) 

