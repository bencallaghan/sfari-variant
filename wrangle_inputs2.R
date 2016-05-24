library(dplyr)
library(ssh.utils)
library(stringr)
library("biomaRt")
library(seqinr)

# GENE <- "SYNGAP1"
# CHROM <- 6
# INPUTDIR <- paste0("/home/bcallaghan/NateDBCopy/inputs/",GENE,"/")
getwd()

## Setup Biomart
listMarts(host="grch37.ensembl.org", path="/biomart/martservice" )
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host="grch37.ensembl.org",path="/biomart/martservice")
listDatasets(mart)
listFilters(mart)
BMattr <- listAttributes(mart)
BMattr[grep("coding",BMattr$name),]
#

## BIOMART - CANONICAL TRANSCRIPT
attributes <- c("ensembl_transcript_id","transcript_start","transcript_end" ,
                "transcript_status","transcript_count", "transcript_biotype", "transcript_source",
                "transcript_version","transcript_length")
BMtranscript <- getBM(attributes = attributes, 
                      filters = c("hgnc_symbol"),values = list(GENE), mart = mart, verbose = FALSE)
BMtranscript
BMtranscript %>% filter(transcript_status == "KNOWN", transcript_biotype == "protein_coding") %>% arrange(desc(transcript_length)) -> BMtranscript.sort
TRANSCRIPT <- BMtranscript.sort$ensembl_transcript_id[1]




## BIOMART - BED
# bed <- getBM(attributes = c("chromosome_name","exon_chrom_start","exon_chrom_end" ,"rank","strand"), 
#              filters = c("hgnc_symbol","refseq_mrna"),values = list("SYNGAP1","NM_006772"), mart = mart, verbose = FALSE)
# colnames(bed) <- c("Chrom","Start","Stop","Exon","Strand")
# bed$Transcript <- TRANSCRIPT
# bed$Gene <- GENE
# bed

BMbed <- getBM(attributes = c("chromosome_name","exon_chrom_start","exon_chrom_end" ,"rank","strand","ensembl_transcript_id"), 
      filters = c("chromosome_name","hgnc_symbol","ensembl_transcript_id"),values = list(CHROM,GENE,TRANSCRIPT), mart = mart, verbose = FALSE)
colnames(BMbed) <- c("Chrom","Start","Stop","Exon","Strand","Transcript")
BMbed$Gene <- GENE
BED <- BMbed
cdna_length <- sum(BED$Stop - BED$Start)
cdna_length



## BIOMART - TRANSCRIPT (cDNA)
BMcdna <- getBM(attributes = c("cdna","refseq_mrna","ensembl_transcript_id"), 
                filters = c("chromosome_name","hgnc_symbol","ensembl_transcript_id"),values = list(CHROM,GENE,TRANSCRIPT), mart = mart, verbose = FALSE)
CDNA <- BMcdna$cdna
cdna <- tolower(CDNA)
nchar(CDNA)
cdna.translated <- translate_cdna(cdna)
cdna.orf <- get_cdna_orf(cdna)
orf.coords <- get_orf_coords(cdna)

## BIOMART - PEPTIDE
BMpeptide <- getBM(attributes = c("peptide","refseq_mrna","ensembl_transcript_id"), 
                   filters = c("hgnc_symbol","ensembl_transcript_id"),values = list(GENE,TRANSCRIPT), mart = mart, verbose = FALSE)
FASTA <- BMpeptide$peptide
FASTA == cdna.translated


## BIOMART - GENOME
BMgene <- getBM(attributes = c("transcript_exon_intron","refseq_mrna","ensembl_transcript_id"), 
                   filters = c("hgnc_symbol","ensembl_transcript_id"),values = list(GENE,TRANSCRIPT), mart = mart, verbose = FALSE)
BMgenecoords <- getBM(attributes = c("start_position","end_position"), 
                      filters = c("hgnc_symbol","ensembl_transcript_id"),values = list(GENE,TRANSCRIPT), mart = mart, verbose = FALSE)

genebases <- data.frame()
# peptideseq <- getBM(attributes = c("peptide","refseq_mrna","ensembl_transcript_id"), 
#                     filters = c("hgnc_symbol","refseq_mrna"),values = list("SYNGAP1","NM_006772"), mart = mart, verbose = FALSE)
# getBM(attributes = c("peptide","refseq_mrna","ensembl_transcript_id","canonical_transcript_stable_id"), filters = c("hgnc_symbol","refseq_mrna"),values = list("SYNGAP1","NM_006772"), mart = mart, verbose = FALSE)


## BED File
# bedpath <- paste0(INPUTDIR,GENE,".bed")
# print(bedpath)
# bedfile <- read.table(bedpath,skip= 1, col.names=c('chr','Start','Stop','transcript','something','strand'))
# exoncounts <- table(gsub("(.+)_exon.+","\\1", bedfile$transcript))
# canonical <- names(exoncounts[which(exoncounts == max(exoncounts))])
# Choose canonical as most exons - if you capture extra (untranscribed) exons it will be filtered out anyway but don't want to lose
# any at this point
# BED <- bedfile[grepl(paste0(canonical,".+"),bedfile$transcript),]
# BED.canon <- read.table('/home/bcallaghan/NateDBCopy/inputs/UCSC_exons_modif_canonical.bed')
# BED.canon %>% filter(V1 == paste0('chr',CHROM) & V4 == GENE) -> BED.canon.gene
# colnames(BED.canon.gene) <- c("Chrom","Start","Stop","Gene","Exon","Transcript","Strand")
# BED <- BED.canon.gene
# 
# bed <- getBM(attributes = c("chromosome_name","exon_chrom_start","exon_chrom_end" ,"rank","strand"), 
#              filters = c("hgnc_symbol","refseq_mrna"),values = list("SYNGAP1","NM_006772"), mart = mart, verbose = FALSE)
# colnames(BED.canon.gene) <- c("Chrom","Start","Stop","Exon","Strand")
# bed$Transcript <- TRANSCRIPT
# bed$Gene <- GENE

##
# 
# ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
# filters = listFilters(ensembl)
# seq = getSequence(id="BRCA1", type="hgnc_symbol", seqType="peptide", mart = ensembl)
# getBM(c("peptide", "hgnc_symbol"))
# 
# 
# listMarts(host="ensembl.org")
# mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
# human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# seq = getSequence(id="BRCA1", type="hgnc_symbol", seqType="peptide", mart = mart)



## FASTA File
# fastapath <- paste0(INPUTDIR,GENE,".fa")
# print(fastapath)
# fastafile <- read.table(fastapath,skip = 1)
# FASTA <- paste0(fastafile$V1,collapse="")

## Exonic CADD file

awkcommand <- scribe_awkcommand(CHROM,BED,GENE,INPUTDIR)
print(awkcommand)

break

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
