library("data.table")
library("ggplot2")

setwd('f:/mmrf_data')

gisticAmpGenesFile <- 'F:/FROM_BRII/mmrf_gistic2_results_2019_0102/mmrf_results/amp_genes.conf_90.txt'


message(paste0('Processing ', basename(gisticAmpGenesFile), '..'))
ampGenes <- data.table::fread(input = gisticAmpGenesFile, stringsAsFactors = FALSE, header = TRUE)

# Removes last column (V51) 
if(colnames(ampGenes)[ncol(ampGenes)] == paste('V', ncol(ampGenes), sep='')){
  ampGenes <- ampGenes[, 1:ncol(ampGenes)-1, with = FALSE]
}

ampGenesStats <- ampGenes[1:3]

cytoband <- colnames(ampGenesStats)[2:length(colnames(ampGenesStats))]

ampGenesStats

amp2 <- transpose(ampGenesStats)

colnames(amp2) <- c("q_value", "residual_q_value", "wide_peak_boundaries")

str(amp2)

amp3 <- amp2[2:nrow(amp2)]

amp3$cytoband = cytoband

head(amp3)

amp3[,c("cytoband", "q_value")]

amp3[,c("chromestr", "peak_range") := tstrsplit(wide_peak_boundaries, ':')]

head(amp3)

amp3[, "chromosome" := substr(chromestr, 4, length(chromestr))]

head(amp3)

amp3[,c("peak_start", "peak_end") := tstrsplit(peak_range, '-')]

amp3[, amp:=1]

head(amp3)

amp4 <- amp3[, .(chromosome, peak_start, peak_end, amp, cytoband, q_value)]

head(amp4)

#we need data only from 4th row
ampGenes <- ampGenes[4:nrow(ampGenes),]
ampGenes$cytoband = gsub(pattern = ' ', replacement = '_', x = ampGenes$cytoband)

data.table::setDF(x = ampGenes)
ampGenes <- suppressWarnings(data.table::melt(ampGenes, id.vars = 'cytoband'))
data.table::setDT(ampGenes)

#Remove empty values
ampGenes <- ampGenes[!value %in% '']

ampGenes <- ampGenes[!grep(pattern = '|', x = ampGenes$value, fixed = TRUE)] #remove genes with ambiguous annotation
ampGenes2 <- ampGenes[, .(variable, value)]

colnames(ampGenes2)

ag3 <- aggregate(ampGenes2, by = list(variable2 = ampGenes2$variable), function(x) paste(x, collapse = ',') )

ag4 <- ag3[, c(1, 3)]

head(ag4)

colnames(ag4) <- c("cytoband", "genes_in_region")

head(ag4)

data.table::setDT(ag4)

amp5 <- amp4[ag4, on="cytoband"]

head(amp5, 30)

colnames(amp5)

amp6 <- amp5[, .(chromosome, peak_start, peak_end, genes_in_region, amp, cytoband, q_value)]

head(amp6,7)

ampOutputFile <- "DT_R/cbio_files/data_gistic2_amp_file.txt"

file.remove(outputFile)
file.create(outputFile)

write.table(amp6, file=ampOutputFile, sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, append = TRUE, na="") 




#---------------------------------------
# DELETIONS FILE FROM GISTIC
#---------------------------------------

gisticDelGenesFile <- 'F:/FROM_BRII/mmrf_gistic2_results_2019_0102/mmrf_results/del_genes.conf_90.txt'

message(paste0('Processing ', basename(gisticDelGenesFile), '..'))
delGenes <- data.table::fread(input = gisticDelGenesFile, stringsAsFactors = FALSE, header = TRUE)

# Removes last column (V51) 
if(colnames(delGenes)[ncol(delGenes)] == paste('V', ncol(delGenes), sep='')){
  delGenes <- delGenes[, 1:ncol(delGenes)-1, with = FALSE]
}

delGenesStats <- delGenes[1:3]

cytoband <- colnames(delGenesStats)[2:length(colnames(delGenesStats))]

delGenesStats

del2 <- transpose(delGenesStats)

colnames(del2) <- c("q_value", "residual_q_value", "wide_peak_boundaries")

str(del2)

del3 <- del2[2:nrow(del2)]

del3$cytoband = cytoband

head(del3)

del3[,c("cytoband", "q_value")]

del3[,c("chromestr", "peak_range") := tstrsplit(wide_peak_boundaries, ':')]

head(del3)

del3[, "chromosome" := substr(chromestr, 4, length(chromestr))]

head(del3)

del3[,c("peak_start", "peak_end") := tstrsplit(peak_range, '-')]

del3[, amp:=0]

head(del3)

del4 <- del3[, .(chromosome, peak_start, peak_end, amp, cytoband, q_value)]

head(del4)

#we need data only from 4th row
delGenes <- delGenes[4:nrow(delGenes),]
delGenes$cytoband = gsub(pattern = ' ', replacement = '_', x = delGenes$cytoband)

data.table::setDF(x = delGenes)
delGenes <- suppressWarnings(data.table::melt(delGenes, id.vars = 'cytoband'))
data.table::setDT(delGenes)

#Remove empty values
delGenes <- delGenes[!value %in% '']

delGenes <- delGenes[!grep(pattern = '|', x = delGenes$value, fixed = TRUE)] #remove genes with ambiguous annotation
delGenes2 <- delGenes[, .(variable, value)]

colnames(delGenes2)

dg3 <- aggregate(delGenes2, by = list(variable2 = delGenes2$variable), function(x) paste(x, collapse = ',') )

dg4 <- dg3[, c(1, 3)]

head(dg4)

colnames(dg4) <- c("cytoband", "genes_in_region")

head(dg4)

data.table::setDT(dg4)

del5 <- del4[dg4, on="cytoband"]

head(del5, 30)

colnames(del5)

del6 <- del5[, .(chromosome, peak_start, peak_end, genes_in_region, amp, cytoband, q_value)]

head(del6,7)

delOutputFile <- "DT_R/cbio_files/data_gistic2_del_file.txt"

file.remove(delOutputFile)
file.create(delOputFile)

write.table(del6, file=delOutputFile, sep = "\t", col.names = TRUE, row.names=FALSE, quote = FALSE, append = TRUE, na="") 




