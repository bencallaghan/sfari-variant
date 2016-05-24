
# FIND starts and stops
sequence <- "aaaaaatga"
findPotentialStartsAndStops <- function(sequence)
{
  # Define a vector with the sequences of potential start and stop codons
  codons            <- c("atg", "taa", "tag", "tga")
  # Find the number of occurrences of each type of potential start or stop codon
  for (i in 1:4)
  {
    codon <- codons[i]
    # Find all occurrences of codon "codon" in sequence "sequence"
    occurrences <- matchPattern(codon, sequence)
    # Find the start positions of all occurrences of "codon" in sequence "sequence"
    codonpositions <- attr(attr(occurrences,"range"),"start")
    # Find the total number of potential start and stop codons in sequence "sequence"
    numoccurrences <- length(codonpositions)
    if (i == 1)
    {
      # Make a copy of vector "codonpositions" called "positions"
      positions <- codonpositions
      # Make a vector "types" containing "numoccurrences" copies of "codon"
      types <- rep(codon, numoccurrences)
    }
    else
    {
      # Add the vector "codonpositions" to the end of vector "positions":
      positions   <- append(positions, codonpositions, after=length(positions))
      # Add the vector "rep(codon, numoccurrences)" to the end of vector "types":
      types       <- append(types, rep(codon, numoccurrences), after=length(types))
    }
  }
  # Sort the vectors "positions" and "types" in order of position along the input sequence:
  indices <- order(positions)
  positions <- positions[indices]
  types <- types[indices]
  # Return a list variable including vectors "positions" and "types":
  mylist <- list(positions,types)
  return(mylist)
}



# Find ORFS
findORFsinSeq <- function(sequence)
{
  require(Biostrings)
  # Make vectors "positions" and "types" containing information on the positions of ATGs in the sequence:
  mylist <- findPotentialStartsAndStops(sequence)
  positions <- mylist[[1]]
  types <- mylist[[2]]
  # Make vectors "orfstarts" and "orfstops" to store the predicted start and stop codons of ORFs
  orfstarts <- numeric()
  orfstops <- numeric()
  # Make a vector "orflengths" to store the lengths of the ORFs
  orflengths <- numeric()
  # Print out the positions of ORFs in the sequence:
  # Find the length of vector "positions"
  numpositions <- length(positions)
  # There must be at least one start codon and one stop codon to have an ORF.
  if (numpositions >= 2)
  {
    for (i in 1:(numpositions-1))
    {
      posi <- positions[i]
      typei <- types[i]
      found <- 0
      while (found == 0)
      {
        for (j in (i+1):numpositions)
        {
          posj  <- positions[j]
          typej <- types[j]
          posdiff <- posj - posi
          posdiffmod3 <- posdiff %% 3
          # Add in the length of the stop codon
          orflength <- posj - posi + 3
          if (typei == "atg" && (typej == "taa" || typej == "tag" || typej == "tga") && posdiffmod3 == 0)
          {
            # Check if we have already used the stop codon at posj+2 in an ORF
            numorfs <- length(orfstops)
            usedstop <- -1
            if (numorfs > 0)
            {
              for (k in 1:numorfs)
              {
                orfstopk <- orfstops[k]
                if (orfstopk == (posj + 2)) { usedstop <- 1 }
              }
            }
            if (usedstop == -1)
            {
              orfstarts <- append(orfstarts, posi, after=length(orfstarts))
              orfstops <- append(orfstops, posj+2, after=length(orfstops)) # Including the stop codon.
              orflengths <- append(orflengths, orflength, after=length(orflengths))
            }
            found <- 1
            break
          }
          if (j == numpositions) { found <- 1 }
        }
      }
    }
  }
  # Sort the final ORFs by start position:
  indices <- order(orfstarts)
  orfstarts <- orfstarts[indices]
  orfstops <- orfstops[indices]
  # Find the lengths of the ORFs that we have
  orflengths <- numeric()
  numorfs <- length(orfstarts)
  for (i in 1:numorfs)
  {
    orfstart <- orfstarts[i]
    orfstop <- orfstops[i]
    orflength <- orfstop - orfstart + 1
    orflengths <- append(orflengths,orflength,after=length(orflengths))
  }
  mylist <- list(orfstarts, orfstops, orflengths)
  return(mylist)
}

plotORFsinSeq <- function(sequence)
{
  # Make vectors "positions" and "types" containing information on the positions of ATGs in the sequence:
  mylist <- findPotentialStartsAndStops(sequence)
  positions <- mylist[[1]]
  types <- mylist[[2]]
  # Make vectors "orfstarts" and "orfstops" to store the predicted start and stop codons of ORFs
  orfstarts <- numeric()
  orfstops <- numeric()
  # Make a vector "orflengths" to store the lengths of the ORFs
  orflengths <- numeric()
  # Print out the positions of ORFs in the sequence:
  numpositions <- length(positions) # Find the length of vector "positions"
  # There must be at least one start codon and one stop codon to have an ORF.
  if (numpositions >= 2)
  {
    for (i in 1:(numpositions-1))
    {
      posi <- positions[i]
      typei <- types[i]
      found <- 0
      while (found == 0)
      {
        for (j in (i+1):numpositions)
        {
          posj <- positions[j]
          typej <- types[j]
          posdiff <- posj - posi
          posdiffmod3 <- posdiff %% 3
          orflength <- posj - posi + 3 # Add in the length of the stop codon
          if (typei == "atg" && (typej == "taa" || typej == "tag" || typej == "tga") && posdiffmod3 == 0)
          {
            # Check if we have already used the stop codon at posj+2 in an ORF
            numorfs <- length(orfstops)
            usedstop <- -1
            if (numorfs > 0)
            {
              for (k in 1:numorfs)
              {
                orfstopk <- orfstops[k]
                if (orfstopk == (posj + 2)) { usedstop <- 1 }
              }
            }
            if (usedstop == -1)
            {
              orfstarts <- append(orfstarts, posi, after=length(orfstarts))
              orfstops <- append(orfstops, posj+2, after=length(orfstops)) # Including the stop codon.
              orflengths <- append(orflengths, orflength, after=length(orflengths))
            }
            found <- 1
            break
          }
          if (j == numpositions) { found <- 1 }
        }
      }
    }
  }
  # Sort the final ORFs by start position:
  indices <- order(orfstarts)
  orfstarts <- orfstarts[indices]
  orfstops <- orfstops[indices]
  # Make a plot showing the positions of ORFs in the input sequence:
  # Draw a line at y=0 from 1 to the length of the sequence:
  x <- c(1,nchar(sequence))
  y <- c(0,0)
  plot(x, y, ylim=c(0,3), type="l", axes=FALSE, xlab="Nucleotide", ylab="Reading frame", main="Predicted ORFs")
  segments(1,1,nchar(sequence),1)
  segments(1,2,nchar(sequence),2)
  # Add the x-axis at y=0:
  axis(1, pos=0)
  # Add the y-axis labels:
  text(0.9,0.5,"+1")
  text(0.9,1.5,"+2")
  text(0.9,2.5,"+3")
  # Make a plot of the ORFs in the sequence:
  numorfs <- length(orfstarts)
  for (i in 1:numorfs)
  {
    orfstart <- orfstarts[i]
    orfstop <- orfstops[i]
    remainder <- (orfstart-1) %% 3
    if    (remainder == 0) # +1 reading frame
    {
      rect(orfstart,0,orfstop,1,col="cyan",border="black")
    }
    else if (remainder == 1)
    {
      rect(orfstart,1,orfstop,2,col="cyan",border="black")
    }
    else if (remainder == 2)
    {
      rect(orfstart,2,orfstop,3,col="cyan",border="black")
    }
  }
}

translate_cdna <- function(cdnaseq){
  forwardORFs <- findORFsinSeq(cdna)
  plotORFsinSeq(cdna)
  cdna.rev <- c2s(rev(comp(s2c(cdna))))
  reverseORFs <- findORFsinSeq(cdna.rev)
  plotORFsinSeq(cdna.rev)
  if(max(forwardORFs[[3]])>max(reverseORFs[[3]])){
    start_trans <- forwardORFs[[1]][which(forwardORFs[[3]] == max(forwardORFs[[3]]))]
    stop_trans <- forwardORFs[[2]][which(forwardORFs[[3]] == max(forwardORFs[[3]]))]
    myorf <- substring(cdna, start_trans, stop_trans)
    cdna.translated <-  c2s(seqinr::translate(s2c(myorf)))                           
  }else{
    stop_trans <- reverseORFs[[1]][which(reverseORFs[[3]] == max(reverseORFs[[3]]))]
    stop_trans <- reverseORFs[[2]][which(reverseORFs[[3]] == max(reverseORFs[[3]]))]
    myorf <- substring(cdna, start_trans, stop_trans)
    cdna.translated <- c2s(seqinr::translate(s2c(myorf)))
    
  }
  return(cdna.translated)
}

get_cdna_orf <- function(cdnaseq){
  forwardORFs <- findORFsinSeq(cdna)
  plotORFsinSeq(cdna)
  cdna.rev <- c2s(rev(comp(s2c(cdna))))
  reverseORFs <- findORFsinSeq(cdna.rev)
  plotORFsinSeq(cdna.rev)
  if(max(forwardORFs[[3]])>max(reverseORFs[[3]])){
    start_trans <- forwardORFs[[1]][which(forwardORFs[[3]] == max(forwardORFs[[3]]))]
    stop_trans <- forwardORFs[[2]][which(forwardORFs[[3]] == max(forwardORFs[[3]]))]
    myorf <- substring(cdna, start_trans, stop_trans)
#     cdna.translated <-  c2s(seqinr::translate(s2c(myorf)))                           
  }else{
    stop_trans <- reverseORFs[[1]][which(reverseORFs[[3]] == max(reverseORFs[[3]]))]
    stop_trans <- reverseORFs[[2]][which(reverseORFs[[3]] == max(reverseORFs[[3]]))]
    myorf <- substring(cdna, start_trans, stop_trans)
#     cdna.translated <- c2s(seqinr::translate(s2c(myorf)))
    
  }
  return(myorf)
}

get_orf_coords <- function(cdnaseq){
  forwardORFs <- findORFsinSeq(cdna)
  plotORFsinSeq(cdna)
  cdna.rev <- c2s(rev(comp(s2c(cdna))))
  reverseORFs <- findORFsinSeq(cdna.rev)
  plotORFsinSeq(cdna.rev)
  if(max(forwardORFs[[3]])>max(reverseORFs[[3]])){
    start_trans <- forwardORFs[[1]][which(forwardORFs[[3]] == max(forwardORFs[[3]]))]
    stop_trans <- forwardORFs[[2]][which(forwardORFs[[3]] == max(forwardORFs[[3]]))]
    myorf <- substring(cdna, start_trans, stop_trans)
    #     cdna.translated <-  c2s(seqinr::translate(s2c(myorf)))                           
  }else{
    stop_trans <- reverseORFs[[1]][which(reverseORFs[[3]] == max(reverseORFs[[3]]))]
    stop_trans <- reverseORFs[[2]][which(reverseORFs[[3]] == max(reverseORFs[[3]]))]
    myorf <- substring(cdna, start_trans, stop_trans)
    #     cdna.translated <- c2s(seqinr::translate(s2c(myorf)))
    
  }
  return(c(start_trans,stop_trans))
}

scribe_awkcommand <- function(CHROM,BED,GENE,INPUTDIR){
  awkcommand <- paste0("awk '$1 == ",CHROM, " && (")
  for (i in 1:nrow(BED)){
    if(i < nrow(BED)){
      print(i)
      exoniccmd <- (paste0("($2 > ",BED$Start[i] - 1, " && $2 < ", BED$Stop[i] + 1 ,") || "))
      awkcommand <- paste0(awkcommand,exoniccmd)
    }else{
      print(i)
      exoniccmd <- (paste0("($2 > ",BED$Start[i] - 1, " && $2 < ", BED$Stop[i] + 1,")) "))
      awkcommand <- paste0(awkcommand,exoniccmd,"{print $0}' /space/bin/annovar/humandb/hg19_cadd.txt > ",
                           INPUTDIR,"cadd",GENE, " 2> ", INPUTDIR, "cadd", GENE, ".err" )
    }
  }
  return(awkcommand)
}

make_anno_file <- function(cdna.orf){
  for(i in 1:length(s2c(cdna.orf))){
    print(s2c(cdna.orf)[i])
  }
  anno.df <- data.frame(matrix(nrow=nchar(cdna.orf),ncol=5))
  colnames(anno.df) <- c("CHR","START","STOP","REF","ALT")
  anno.df$REF <- s2c(cdna.orf)
}
