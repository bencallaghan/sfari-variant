GENE <- "SYNGAP1"
CHROM <- 6
# TRANSCRIPT <-  "NM_006772"

INPUTDIR <- paste0("/home/bcallaghan/NateDBCopy/inputs/",GENE,"/")
getwd()

# Source scripts
source("functions.R")
source("wrangle_inputs2.R")
source("GENE_prioritised_variants2.R")
