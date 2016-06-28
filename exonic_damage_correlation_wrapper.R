setwd("/home/bcallaghan/NateDBCopy/20kcorr/")
predictProteinMap <- read.table("fastamap")

for(i in nrow(predictProteinMap)){
  GENE = predictProteinMap$V2[i]
  source("/home/bcallaghan/NateDBCopy/exonic_damage_correlation.R")
}