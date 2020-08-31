
setwd("f:/mmrf_data/DT_R")

# Run Scripts to generate cbio_files

source("DT_CREATE_CLINICAL_PATIENT_FILE.r")

setwd("f:/mmrf_data/DT_R")

source("DT_CREATE_CLINICAL_SAMPLE_FILE.r")

setwd("f:/mmrf_data/DT_R")

# CREATES MUTATION DATA FILE AND case_lists/cases_sequenced.txt
source("DT_CREATE_MUTATION_FILE.r")

setwd("f:/mmrf_data/DT_R")

# CREATES EXPRESSION DATA FILE AND case_lists/cases_rna_seq_mrna.txt
source("DT_CREATE_EXPRESSION.r")

setwd("f:/mmrf_data/DT_R")

# CREATES Discrete CNA DATA FILE AND case_lists/cases_cna.txt
source("DT_CREATE_DISCRETE_COPY_NUMBER_DATA.r")

setwd("f:/mmrf_data/DT_R")

tar("cbio_files.tgz", "cbio_files", compression="gzip")

