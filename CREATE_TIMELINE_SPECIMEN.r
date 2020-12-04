library(purrr)
library(dplyr)
library(tidyr)

setwd('c:/Users/abhmalat/OneDrive - Indiana University/MMRF_CoMMpass_IA16a')

specimen <- read.csv('clinical_flat_files/CoMMpass_IA16_FlatFiles/MMRF_CoMMpass_IA16_PER_PATIENT_VISIT.csv',
                       sep=",", header=TRUE) %>%  dplyr::select(c("PUBLIC_ID", "VISITDY", "Specimen_ID", "Sample_Type"))

specimen <- specimen[!(specimen$Specimen_ID == ""),]
specimen$EVENT_TYPE <- "SPECIMEN"
specimen$STOP_DATE <- ""

names(specimen) <- c("PATIENT_ID", "START_DATE", "SPECIMEN_REFERENCE_NUMBER", "Sample Type", "EVENT_TYPE", "STOP_DATE")

outputFile <- "data_timeline_specimen.txt"
write.table(specimen %>% dplyr::select(c("PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "SPECIMEN_REFERENCE_NUMBER", "Sample Type")),
                                       outputFile, sep='\t', col.names = TRUE, row.names=FALSE, quote = FALSE, append = FALSE)

head(per_visit)

bone_assessment <- per_visit[, c("PUBLIC_ID", "BA_WASABONEASSES", "BA_DAYOFASSESSM", "BA_TYPEOFBONEASS", "BA_SPECIFY", "BA_IMAGEINTERPRE",
               "BA_LYTICLESIONS", "BA_OFLYTICLESION", "BA_OSTEOPENIAOST", "BA_PATHOLOGICFRA", "BA_MULTIPLEDISSE",
               "BA_SKULL", "BA_SPINE", "BA_CERVICAL", "BA_THORACIC", "BA_LUMBAR", "BA_SACRAL", "BA_PELVIS", "BA_RIGHT",
               "BA_LEFT", "BA_FEMUR", "BA_RIGHT2", "BA_LEFT2", "BA_HUMERUS", "BA_RIGHT3", "BA_LEFT3", "BA_RIBS", "BA_RIGHT4",
               "BA_NUMBERAFFECTE", "BA_LEFT4", "BA_NUMBERAFFECTE2", "BA_SPECIFY2")]

ba_true <- dplyr::filter(bone_assessment, BA_WASABONEASSES == "Yes")

ba_true$BA_TYPEOFBONEASS <- ifelse(ba_true$BA_TYPEOFBONEASS %in% "Other", ba_true$BA_SPECIFY, ba_true$BA_TYPEOFBONEASS)
ba_true <- dplyr::select (ba_true,-c(BA_SPECIFY, BA_WASABONEASSES, ))
head(ba_true)

ba_true$STOP_DATE <- ""
ba_true$EVENT_TYPE <- "LAB_TEST"
ba_true$TEST <- "Bone Assessment"

colnames(ba_true)
names(ba_true) <- c("PATIENT_ID", "START_DATE", "Type of bone assessment", "Image Interpretation", "Lytic Lesions",
               "Number of Lytic Lesions", "Osteopenia / Osteoporosis", "Pathologic fracture", "Multiple disseminated lesions",
               "Skull", "Spine", "Cervical", "Thoracic", "Lumbar", "Sacral", "Pelvis", "Right pelvis", "Left pelvis",
               "Femur", "Right femur", "Left femur", "Humerus", "Right humerus", "Left humerus", "Ribs", "Right ribs",
               "Number of ribs affected (right side)", "Left ribs", "Number of ribs affected (left side)",
               "Other myeloma-related bone lesion", "STOP_DATE", "EVENT_TYPE", "TEST")

summary(bone_assessment)

outputFile <- "data_timeline_lab_bone_assessment.txt"

print('Completed')