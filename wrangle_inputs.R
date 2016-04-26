library(dplyr)
library(ssh.utils)
library(stringr)

GENE <- "SYNGAP1"
CHROM <- 6
INPUTDIR <- paste0("/home/bcallaghan/NateDBCopy/inputs/",GENE,"/")
getwd()

#

## BED File
bedpath <- paste0(INPUTDIR,GENE,".bed")
print(bedpath)
bedfile <- read.table(bedpath,skip= 1, col.names=c('chr','Start','Stop','transcript','something','strand'))
exoncounts <- table(gsub("(.+)_exon.+","\\1", bedfile$transcript))
canonical <- names(exoncounts[which(exoncounts == max(exoncounts))])
# Choose canonical as most exons - if you capture extra (untranscribed) exons it will be filtered out anyway but don't want to lose
# any at this point
BED <- bedfile[grepl(paste0(canonical,".+"),bedfile$transcript),]
##

## FASTA File
fastapath <- paste0(INPUTDIR,GENE,".fa")
print(fastapath)
fastafile <- read.table(fastapath,skip = 1)
FASTA <- paste0(fastafile$V1,collapse="")

## Exonic CADD file
awkcommand <- paste0("awk '$1 == ",CHROM, " && (")
for (i in 1:nrow(BED)){
  if(i < nrow(BED)){
  print(i)
  exoniccmd <- (paste0("($2 > ",BED$Start[i], " && $2 < ", BED$Stop[i],") || "))
  awkcommand <- paste0(awkcommand,exoniccmd)
  }else{
    print(i)
    exoniccmd <- (paste0("($2 > ",BED$Start[i], " && $2 < ", BED$Stop[i],")) "))
    awkcommand <- paste0(awkcommand,exoniccmd,"{print $0}' /space/bin/annovar/humandb/hg19_cadd.txt > ", INPUTDIR,"cadd",GENE )
  }
}
print(awkcommand)
touchcmd <- paste0("touch ", INPUTDIR, "cadd",GENE)
awkcmd <- paste0(touchcmd, " ; ", awkcommand)
# cmd.out <- run.remote(cmd=awkcmd , remote= "apu")

## Run annovar on exonic CADD
caddpath <- paste0(INPUTDIR,"cadd",GENE)
annocmd <- paste0("perl /space/bin/annovar/table_annovar.pl ", caddpath, " /space/bin/annovar/humandb/ -buildver hg19 -out ", INPUTDIR,GENE, "_anno -otherinfo -remove -protocol refGene,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp135,ljb_all,exac03,cadd -operation g,r,f,f,f,f,f,f -nastring . -csvout")
# cmd.out <- run.remote(cmd=annocmd , remote= "apu")

## File checks:
check_inputs <- function(){
  miss <- NULL
  if(file.exists(fastapath) == FALSE){
    miss <- "fasta"
  } 
  if(file.exists(bedpath) == FALSE){
    miss <- c(miss,"bed")
  } 
  if(file.exists(caddpath) == FALSE){
    miss <- c(miss,"cadd")
  } 
  if(file.exists(paste0(INPUTDIR,GENE,"_anno.hg19_multianno.csv")) == FALSE){
    miss <- c(miss,"caddanno")
  } 
  if(is.null(miss)){
    return("all good")
  }else{
  return(paste(miss, "files missing, fix before prioritisation."))
}
}

check_inputs()
