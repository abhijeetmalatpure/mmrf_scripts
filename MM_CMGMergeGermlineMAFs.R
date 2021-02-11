library(dplyr)
library(stringr)
library(data.table)
library(maftools)

# Old version that reads them all and then tries to merge
# combine_maf_files <- function (filenames) {
#   maf_files <- list()
#   failed_maf <- list()
#   x <- 1
#   for (i in seq_along(filenames)) {
#     print(paste(i,"-", "Reading", filenames[i]))
#     maf <- read_maf_file(filenames[i])#, vc_nonSyn = vcNames)
#     if (!is.na(maf)) {
#       maf_files[x] <- maf
#       x <- x + 1 }
#     else {
#       failed_maf <- append(failed_maf, filenames[i])
#     }
#   }
#   print("***********************************************************************")
#   return(list("merged"= merge_mafs(maf_files), "failed" = failed_maf))
# }

# New version that sequentially reads and merges as it proceeds
combine_maf_files <- function (filenames) {
  for (i in seq_along(filenames)) {
    print(paste(i,"-", "Reading", filenames[i]))
    maf <- read_maf_file(filenames[i])#, vc_nonSyn = vcNames)
    if (!is.na(maf)) {
      if (i == 1) {
        final_maf <- maf
        print("Reading next file to merge into first..")
      }
      else {
        print(paste(filenames[i], "read correctly. Merging into the previous: "))
        final_maf <- merge_maf_files(final_maf, maf)
      }
    }
    else {
      failed_maf <- append(failed_maf, filenames[i])
    }
  }
  return(list("merged"= final_maf, "failed" = failed_maf))
}


# Read MAF file using read.maf. Return NaN if it errors out
read_maf_file <- function(x) {
  maf <- data.table::fread(x)
  vcNames <- unique(maf$Variant_Classification)
  tryCatch(read.maf(x, vc_nonSyn = vcNames),
           error = function(x) {print(paste(x, ": Error reading file. Skipping")); NaN})
}


# Merge two MAF files into one. Return first if merging error occurs
merge_maf_files <- function(maf_merged, maf_b) {
  tryCatch(merge_mafs(c(maf_merged, maf_b)),
           error = function(x) {print(paste(": Error merging this file:", x, "Skipping")); maf_merged})
}

mafdir <- "/N/slate/abhmalat/cmg_mm/vcf2mafConversion/maf/germline"
setwd(mafdir)

germlinemafs <- list.files(mafdir, pattern = "*\\.germline.maf$", recursive = FALSE)
germlinemafs
# Call combine_maf_files
processed <- combine_maf_files(germlinemafs)

# Save merged MAF file and failed list of MAFs
merged <- processed$merged
failed <- processed$failed

merged@data$Mutation_Status <- as.character(merged@data$Mutation_Status)
merged@data$Mutation_Status <- "Germline" # Or "Germline" for germline


# Write merged MAF file
write.table(as.data.frame(merged@data), "mm_cmg_germline.maf", col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE, append = FALSE, na = "")

# Save failed maf file list
write.table(as.data.frame(failed), "mm_cmg_germline_failed.txt", col.names = FALSE, row.names = FALSE, sep = "\n", quote = FALSE, append = FALSE, na = "")

mutationCL <- ("cases_sequenced.txt")

f <- file(mutationCL)
writeLines(c(
  "cancer_study_identifier: PHI_MM_2020",
  "stable_id: PHI_MM_2020_sequenced",
  "case_list_name: Sequenced Tumors",
  "case_list_description: All Sequenced Tumors",
  paste("case_list_ids: ", paste(unique(merged@data$Tumor_Sample_Barcode), collapse = '\t'))
  ), f
)
print('MM Germline Mutation case list file completed.')
close(f)
