library(purrr)
library(dplyr)

setwd('c:/Users/abhmalat/OneDrive - Indiana University/MMRF_CoMMpass_IA16a')

#bone_cols <- c("PUBLIC_ID", "BA_WASABONEASSES", "BA_DAYOFASSESSM", "BA_TYPEOFBONEASS", "BA_SPECIFY", "BA_IMAGEINTERPRE",
#               "BA_LYTICLESIONS", "BA_OFLYTICLESION", "BA_OSTEOPENIAOST", "BA_PATHOLOGICFRA", "BA_MULTIPLEDISSE",
#               "BA_SKULL", "BA_SPINE", "BA_CERVICAL", "BA_THORACIC", "BA_LUMBAR", "BA_SACRAL", "BA_PELVIS", "BA_RIGHT",
#               "BA_LEFT", "BA_FEMUR", "BA_RIGHT2", "BA_LEFT2", "BA_HUMERUS", "BA_RIGHT3", "BA_LEFT3", "BA_RIBS", "BA_RIGHT4",
#               "BA_NUMBERAFFECTE", "BA_LEFT4", "BA_NUMBERAFFECTE2", "BA_OTHER", "BA_SPECIFY2")

per_visit <- read.csv('clinical_flat_files/CoMMpass_IA16_FlatFiles/MMRF_CoMMpass_IA16_PER_PATIENT_VISIT.csv', sep=",", header=TRUE)

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

#write.table(select(ba_true, "PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "TEST",
#               "Type of bone assessment", "Image Interpretation", "Lytic Lesions",
#               "Number of Lytic Lesions", "Osteopenia / Osteoporosis", "Pathologic fracture", "Multiple disseminated lesions",
#               "Skull", "Spine", "Cervical", "Thoracic", "Lumbar", "Sacral", "Pelvis", "Right pelvis", "Left pelvis",
#               "Femur", "Right femur", "Left femur", "Humerus", "Right humerus", "Left humerus", "Ribs", "Right ribs",
#               "Number of ribs affected (right side)", "Left ribs", "Number of ribs affected (left side)",
#               "Other myeloma-related bone lesion"),
#            outputFile, sep='\t', col.names = TRUE, row.names=FALSE, quote = FALSE, append = FALSE, na="NA")

mprotein <- per_visit[, c("PUBLIC_ID", "VISIT", "AT_INCREASEOF25F", "AT_SERUMMCOMPONE", "AT_URINEMCOMPONE","AT_ONLYINPATIENT","AT_ONLYINPATIENT2")] %>% dplyr::filter(AT_INCREASEOF25F == "Checked")
mprotein$STOP_DATE <- ""
mprotein$EVENT_TYPE <- "LAB_TEST"
mprotein$TEST <- ""

#mprotein <- dplyr::select (mprotein,-AT_INCREASEOF25F)

#Serum m protein
serum <- select(mprotein, c("PUBLIC_ID", "VISIT", "STOP_DATE", "EVENT_TYPE", "TEST", "AT_SERUMMCOMPONE")) %>% dplyr::filter(AT_SERUMMCOMPONE == "Checked")
serum$TEST <- "Serum M-Component"
names(serum) <- c("PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "TEST", "Serum M-component (absolute increase must be >= 0.5g/dL)")
outputFile <- "data_timeline_lab_m_serum.txt"
write.table(serum, outputFile, sep='\t', col.names = TRUE, row.names=FALSE, quote = FALSE, append = FALSE, na="NA")
metaFile <- file("meta_timeline_lab_m_serum.txt")
writeLines(c(
  "cancer_study_identifier: mmrf_2020",
  "genetic_alteration_type: CLINICAL",
  "datatype: TIMELINE",
  paste("data_filename: ", outputFile)
), metaFile
)
close(metaFile)

#Urine m protein
urine <- select(mprotein, c("PUBLIC_ID", "VISIT", "STOP_DATE", "EVENT_TYPE", "TEST", "AT_URINEMCOMPONE")) %>% dplyr::filter(AT_URINEMCOMPONE == "Checked")
urine$TEST <- "Urine M-Component"
names(urine) <- c("PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "TEST", "Serum M-component (absolute increase must be >= 0.5g/dL)")
outputFile <- "data_timeline_lab_m_urine.txt"
write.table(urine, outputFile, sep='\t', col.names = TRUE, row.names=FALSE, quote = FALSE, append = FALSE, na="NA")
metaFile <- file("meta_timeline_lab_m_urine.txt")
writeLines(c(
  "cancer_study_identifier: mmrf_2020",
  "genetic_alteration_type: CLINICAL",
  "datatype: TIMELINE",
  paste("data_filename: ", outputFile)
), metaFile
)
close(metaFile)

#FLC
flc_bmpcp <- select(mprotein, c("PUBLIC_ID", "VISIT", "STOP_DATE", "EVENT_TYPE", "TEST", "AT_ONLYINPATIENT", "AT_ONLYINPATIENT2")) %>% dplyr::filter(AT_ONLYINPATIENT == "Checked" | AT_ONLYINPATIENT2 == "Checked")
flc_bmpcp$TEST <- "FLC and Bone Marrow Plasma %"
names(flc_bmpcp) <- c("PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "TEST",
                  "Only in patients without measurable serum and M-Protein, the difference between involved and uninvolved FLC levels (absolute increase must be > 10mg/dL)",
                  "Only in patients without measurable serum and urine M-protein and without measurable disease by FLC levels, bone marrow plasma cell percentage (absolute % must be >= 10%)"
)
outputFile <- "data_timeline_lab_flc.txt"
write.table(flc_bmpcp, outputFile, sep='\t', col.names = TRUE, row.names=FALSE, quote = FALSE, append = FALSE, na="NA")
metaFile <- file("meta_timeline_lab_flc.txt")
writeLines(c(
  "cancer_study_identifier: mmrf_2020",
  "genetic_alteration_type: CLINICAL",
  "datatype: TIMELINE",
  paste("data_filename: ", outputFile)
), metaFile
)
close(metaFile)

print('Completed')