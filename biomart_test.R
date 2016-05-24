source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)

listMarts(host="grch37.ensembl.org", path="/biomart/martservice" )

mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host="grch37.ensembl.org",path="/biomart/martservice")
listDatasets(mart)
listFilters(mart)
listAttributes(mart)

affyids=c("202763_at","209310_s_at","207500_at")
getBM(attributes=c('affy_hg_u133_plus_2'  ,  'hgnc_id'), 
      filters =  'affy_hg_u133_plus_2', values = affyids, mart = mart)

entrez=c("673","7157","837")

getsequence <- function(id,seqType,type,mart,verbose = TRUE){
  sequence = getBM(c(seqType, type), filters = type, 
                 values = id, mart = mart, verbose = verbose)
}

getsequence(id = entrez, type="entrezgene",seqType="peptide", mart=mart)

sequence <- getBM(attributes = c("peptide","entrezgene"), filters = "entrezgene",values = c("673"), mart = mart, verbose = FALSE)


bed <- getBM(attributes = c("chromosome_name","exon_chrom_start","exon_chrom_end" ,"rank","strand"), 
             filters = c("hgnc_symbol","refseq_mrna"),values = list("SYNGAP1","NM_006772"), mart = mart, verbose = FALSE)
colnames(BED.canon.gene) <- c("Chrom","Start","Stop","Exon","Strand")
bed$Transcript <- TRANSCRIPT
bed$Gene <- GENE



peptideseq <- getBM(attributes = c("peptide","refseq_mrna","ensembl_transcript_id"), filters = "hgnc_symbol",values = c("SYNGAP1"), mart = mart, verbose = FALSE)
peptideseq <- getBM(attributes = c("peptide","refseq_mrna","ensembl_transcript_id"), filters = c("hgnc_symbol","refseq_mrna"),values = list("SYNGAP1","NM_006772"), mart = mart, verbose = FALSE)

cdnaseq <- getBM(attributes = c("cdna","refseq_mrna","ensembl_transcript_id"), filters = "hgnc_symbol",values = c("SYNGAP1"), mart = mart, verbose = FALSE)
cdnaseq <- getBM(attributes = c("cdna","refseq_mrna","ensembl_transcript_id"), filters = c("hgnc_symbol","refseq_mrna"),values = list("SYNGAP1","NM_006772"), mart = mart, verbose = FALSE)



#If proxy doesn't work...?
# Sys.putenv("http\_proxy" = "http://my.proxy.org:9999")

