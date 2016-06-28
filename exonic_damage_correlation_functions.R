bioGetTranscript<- function(GENE){
  attributes <- c("refseq_mrna","ensembl_transcript_id","transcript_start","transcript_end" ,
                  "transcript_status","transcript_count", "transcript_biotype", "transcript_source",
                  "transcript_version","transcript_length")
  BMtranscript <- getBM(attributes = attributes, 
                        filters = c("hgnc_symbol"),values = list(GENE), mart = mart, verbose = FALSE)
  BMtranscript %>% filter(transcript_status == "KNOWN", transcript_biotype == "protein_coding",refseq_mrna !="") %>% 
    arrange(desc(transcript_length)) -> BMtranscript.sort
}
bioGetBED <- function(GENE){
  
}

#_---3-3-3-3-3-3-3-3--3-3-3-3-3-3-3-3-3-3--3-3-3
retrieve_exonic_vcf <- function(GENE){
  # Take in gene, match to Biomart
  # Find canonical transcript (transcript of largest size)
 
  BMtranscript.sort <- bioGetTranscript(GENE)
  ensembl_ID <- BMtranscript.sort$ensembl_transcript_id[1]
  refseq_ID <- BMtranscript.sort$refseq_mrna[1]
  if(nrow(BMtranscript.sort) == 0){
    return(FALSE) # If not in BM, exit early and return blank
  }
  

  # Find Chromosome
  # Retrieve exonic BED for canonical transcript
  BED <- getBM(attributes = c("chromosome_name","exon_chrom_start","exon_chrom_end" ,"rank","strand","ensembl_transcript_id"), 
               filters = c("hgnc_symbol","ensembl_transcript_id"),values = list(GENE,ensembl_ID), mart = mart, verbose = FALSE)
  CHROM <- BED$chromosome_name[1]
  colnames(BED) <- c("Chrom", "Start", "Stop", "rank", "strand", "ensembl_transcript_id")
  BED$Gene <- GENE
  
  cdna_length <- sum(BED$Stop - BED$Start)
  
  
  # Get exonic DNA sequence
  BMgene <- getBM(attributes = c("transcript_exon_intron","refseq_mrna","ensembl_transcript_id"), 
                  filters = c("hgnc_symbol","ensembl_transcript_id"),values = list(GENE,ensembl_ID), mart = mart, verbose = FALSE)
  BMgenecoords <- getBM(attributes = c("start_position","end_position","genomic_coding_start","genomic_coding_end", "ensembl_transcript_id"), 
                        filters = c("hgnc_symbol","ensembl_transcript_id"),values = list(GENE,ensembl_ID), mart = mart, verbose = FALSE)
  
#   cdna.orf
  genomic_dna <- BMgene$transcript_exon_intron
  # Build VCF
  vcf <- genomic_dna_for_annovar(BMgene$transcript_exon_intron,BED)
  # Get CDNA sequence
  return(vcf)
}


merge_multianno_and_snap2 <- function(multiannodf, snap2df){
  # Merges annovar output and snap2 results
  # Filters for exonic variants in canonical transcript (by filterGenomicVariants())
  # Sanity checks for length and mapping of amino acids
  print("Merging annovar and snap2 results...")
  head(multiannodf)
  head(snap2df)
  
  multiannodf <- filterGenomicVariants(multiannodf,GENE,TRANSCRIPT)
  snap2df$V1 <- levels(snap2df$V1)[snap2df$V1]
  
  cat(nrow(snap2df)/19, "Predict Protein amino acids\n")
  cat(snap2df$V1[1], "...",snap2df$V1[nrow(snap2df)],"\n")
  cat(nrow(multiannodf)/9, "Isoform amino acids\n")
  cat(multiannodf$aachange[1], "...",multiannodf$aachange[nrow(multiannodf)],"\n")
  cat(sum(multiannodf$aachange %in% snap2df$V1)/nrow(multiannodf)*100,"% Amino acid Alignment\n")# Predict Protein doesn't score stopgains/losses
  
  multiannodf$CADD.phred <- phredScale2(multiannodf$cadd)
  snap2df$V3 <- gsub("[A-Z]([0-9]+)[A-Z]","\\1",snap2df$V1)
  
  colnames(snap2df) <- c("snap2aachange","snap2","snap2aapos")
  mergedDF <- merge(multiannodf, snap2df, by.x="aachange", by.y="snap2aachange",all.x=FALSE) 
  mergedDF %>% arrange((Start)) -> mergedDF
  

    
  return(mergedDF)
}

plot_correlation_stuff <- function(mergedDF, mode = "correlation"){
  corr <- function(mergedDF){
    outpath <- paste0(INPUTDIR, "outputs/", GENE, "corr.png")
    p1 <- ggplot(mergedDF,aes(x= CADD.phred, y = snap2)) + geom_point()
    ggsave(outpath,p1)
  }
  snap2 <- function(mergedDF){
    outpath <- paste0(INPUTDIR, "outputs/", GENE, "snap2.png")
    p1 <- ggplot(mergedDF,aes(x= CADD.phred, y = snap2)) + geom_point()
    ggsave(outpath,p1)
  }
  switch(mode, correlation = corr())
  outpath <- paste0(INPUTDIR, "outputs/", GENE, "_CScor.png")
  p1 <- ggplot(multianno.2,aes(x= CADD.phred, y = snap2)) + geom_point()
  ggsave(pngoutputpath, p1)
}