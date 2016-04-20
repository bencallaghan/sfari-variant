library(dplyr)
library(stringr)

GENE <- "SYNGAP1"
getwd()

## BED File
bedpath <- paste0("./inputs/",GENE,"/",GENE,".bed")
print(bedpath)
bedfile <- read.table(bedpath,skip= 1, col.names=c('chr','start','stop','transcript','something','strand'))
exoncounts <- table(gsub("(.+)_exon.+","\\1", bedfile$transcript))
exoncounts[which(max(exoncounts))]

# Choose canonical as most exons - if you capture extra (untranscribed) exons it will be filtered out anyway but don't want to lose
# any at this point

