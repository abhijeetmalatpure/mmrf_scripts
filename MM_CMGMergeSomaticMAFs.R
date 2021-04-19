library(dplyr)
library(stringr)
library(data.table)
library(maftools)

# New version that sequentially reads and merges as it proceeds
# combine_maf_files <- function (filenames) {
#   failed_maf <- c()
#   for (i in seq_along(filenames)) {
#     print(paste(i,"-", "Reading", filenames[i]))
#     maf <- read_maf_file(filenames[i])#, vc_nonSyn = vcNames)
#     if (!is.na(maf)) {
#       if (!exists("final_maf")) {
#         final_maf <- maf
#         print("Reading next file to merge into first..")
#       }
#       else {
#         print(paste(filenames[i], "read correctly. Merging into the previous: "))
#         merged_maf <- merge_maf_files(final_maf, maf)
#         if(!is.na(merged_maf)) {
#           final_maf <- merged_maf
#         }
#         else {
#           print(paste(filenames[i], "could not be merged."))
#           failed_maf <- append(failed_maf, paste0('MERGE_ERR_',filenames[i]))
#         }
#         rm(maf)
#         # Write results to file after every 10 MAF file merges
#         # if(i%%10 == 0) {
#         #   #intermediate_file <- paste0("intermediate_", i, ".maf")
#         #   intermediate_file <- paste0("final_", i, ".maf")
#         #   print(paste("Writing intermediate results to file: ", intermediate_file))
#         #   write.table(as.data.frame(final_maf@data), intermediate_file, col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE, append = FALSE, na = "")
#         #   write.table(as.data.frame(failed_maf), "mm_cmg_germline_failed.txt", col.names = FALSE, row.names = FALSE, sep = "\n", quote = FALSE, append = FALSE, na = "")
#         #   rm(final_maf)
#         # }
#       }
#     }
#     else {
#       failed_maf <- append(failed_maf, paste0('READ_ERR_',filenames[i]))
#     }
#     # for (itm in ls()) {
#     #   print(paste(itm, ":", format(object.size(get(itm)), units = "auto")))
#     # }
#     gc(verbose = TRUE, full = TRUE)
#   }
#   return(list("merged"= final_maf, "failed" = failed_maf))
# }


# Read MAF file using read.maf. Return NaN if it errors out
# read_maf_file <- function(x) {
#   maf <- data.table::fread(x)
#   vcNames <- unique(maf$Variant_Classification)
#   tryCatch(read.maf(x, vc_nonSyn = vcNames),
#            error = function(x) {print(paste(x, ": Error reading file. Skipping")); NA})
# }


# Merge two MAF files into one. Return first if merging error occurs
# merge_maf_files <- function(maf_merged, maf_b) {
#   tryCatch(merge_mafs(c(maf_merged, maf_b)),
#            error = function(x) {print(paste(": Error merging this file:", x, "Skipping")); remove(maf_b); gc(); NA})
# }


combine_maf_files <- function(filenames) {
  final_maf <- c()
  for(f in filenames) {
      maf <- read.csv(f, sep="\t", header = TRUE, comment.char = "#")
    if(length(final_maf)==0) {
      final_maf <- maf
    }
    else {
      final_maf <- rbindlist(list(final_maf, maf), fill = TRUE)
    }
  }
  return(final_maf)
}

mafdir <- "/N/u/abhmalat/Carbonate/slate/cmg_mm/maf/somatic/remerge"
setwd(mafdir)

#germlinemafs <- list.files(mafdir, pattern = "*\\.germline.maf$", recursive = FALSE)
pattrn <- ""#"indels"#"snvs"# "somaticSV"
somaticmafs <- file.path(list.files(mafdir, recursive = FALSE))
somaticmafs
# Call combine_maf_files
processed <- combine_maf_files(somaticmafs)

gc()
# Save merged MAF file and failed list of MAFs
#merged <- processed$merged
#failed <- processed$failed

#merged@data$Mutation_Status <- as.character(merged@data$Mutation_Status)
processed$Mutation_Status <- as.character(processed$Mutation_Status)
processed$Mutation_Status <- "Somatic" # Or "Germline" for germline
processed$NCBI_Build <- "GRCh38"

# Write merged MAF file
write.table(as.data.frame(processed), paste0("mm_cmg_somatic", pattrn, ".maf"), col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE, append = FALSE, na = "")

# Save failed maf file list
#write.table(as.data.frame(failed), paste0("mm_cmg_somatic_", pattrn, "_failed.txt"), col.names = FALSE, row.names = FALSE, sep = "\n", quote = FALSE, append = FALSE, na = "")
#
# mutationCL <- ("cases_sequenced.txt")
#
# f <- file(mutationCL)
# writeLines(c(
#   "cancer_study_identifier: PHI_MM_2020",
#   "stable_id: PHI_MM_2020_sv",
#   "case_list_name: Sequenced Tumors",
#   "case_list_description: All Sequenced Tumors",
#   paste("case_list_ids: ", paste(unique(merged@data$Tumor_Sample_Barcode), collapse = '\t'))
#   ), f
# )
print('MM Germline Mutation case list file completed.')
#close(f)
