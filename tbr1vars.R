library(dplyr)
library(xtable)
library(gridExtra)
library(ggplot2)

setwd("/home/bcallaghan/NateDBCopy")
# png(paste0("out/",sample,".png"), width = 1500, height= (25 * (nrow(p)+1)))
# grid.arrange(tableGrob(p))
# dev.off()
natevars <- read.csv(file="outputs/TBR1/nateANNO.hg19_multianno.csv")
head(natevars)
natevars %>% select(c(1,2,3,4,5,6,7,9,10,26,27)) -> p 

png(filename="outputs/TBR1/natvars.png",width=1500, height = 26 * nrow(p))
grid.arrange(tableGrob(p))
dev.off()

exacvars <- read.csv("outputs/TBR1/TBR1_exacANNO.csv")

exacvars %>% select(c(1:7,9,10,26,27))->p

png("outputs/TBR1/exacvars.png", width = 1500, height = 26 * nrow(p))
grid.arrange(tableGrob(p))
dev.off()

#### 
tbr1 <- read.csv("/home/bcallaghan/awkout_anno.hg19_multianno.csv")
tbr1PP <- read.table("/home/bcallaghan/NateDBCopy/tbr1PP.tsv")

# Fix anno version up
x <- tbr1$Otherinfo
y <- strsplit(as.character(x), split = '\t')
y <- unlist(y)
y <- as.numeric(y[seq(2,length(y),2)])
tbr1$CADD.phred <- as.numeric(as.character(y))
# tbr1 %>% filter(Func.refGene == "exonic") -> tbr1 #Just interested in the exonic...for now
tbr1 %>% mutate(aapos = as.numeric(gsub("TBR1:NM_[0-9]+:exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.[A-Z]([0-9]+)[A-Z]","\\1",tbr1$AAChange.refGene))) -> tbr1


ggplot(tbr1, aes(x= aapos, y = CADD.phred)) + geom_point()
tbr1$aachange <- gsub("TBR1:NM_[0-9]+:exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.","",tbr1$AAChange.refGene,perl=TRUE)

unique(gsub("[A-Z]([0-9]+)[A-Z]","\\1",tbr1$aachange)) #missing exon....

# What's PredictProtein score look like in the context of a gene?
tbr1PP %>% mutate(aapos = as.numeric(gsub('[A-Z]([0-9]+)[A-Z]','\\1',V1))) -> tbr1PP
ggplot(tbr1PP, aes(x = aapos, y=V2)) + geom_point()

# Merge them
tbr1PP$V1 <- as.character(tbr1PP$V1)
tbr1$aachange
tbr1.mr <- merge(tbr1,tbr1PP, by.x='aachange', by.y = 'V1', all=FALSE)
# names(tbr1.mr)[31] <- "PredictProtein"
tbr1.mr$CADD.phred <- as.numeric(tbr1.mr$CADD.phred)
tbr1.mr$PredictProtein <- as.numeric(tbr1.mr$V2)
tbr1.mr$exac03 <- as.numeric(as.character(tbr1.mr$exac03))
tbr1.mr$LJB_MutationTaster <- as.numeric(as.character(tbr1.mr$LJB_MutationTaster))
tbr1.mr$exac03[is.na(tbr1.mr$exac03)] <- 0
# tbr1.mr %>% select(c(1:11,25:31)) %>% arrange(Start) -> tbr1.mr
tbr1.mr %>% arrange(Start) -> tbr1.mr

## Filter Positive Controls
tbr1.mr %>% mutate(PC = ifelse(PredictProtein > 75 & CADD.phred > 25 & (exac03 == 0 | is.na(exac03)) ,yes="PositiveControl", no = FALSE)) -> tbr1.mr
# tbr1.mr %>% mutate(NC = ifelse(PredictProtein > 50 & CADD.phred > 25 & exac03 > 0 ,yes=TRUE, no = FALSE)) -> tbr1.mr

# ggplot(tbr1.mr, aes(x = aapos.x, y = PredictProtein)) + geom_point(aes(colour=PC))# Where is that exon?????
# Graph the Positive Controls
p1 <- ggplot(tbr1.mr, aes(x = aapos.x, y = PredictProtein)) + geom_point(aes(colour=PC)) + xlab("Amino Acid Coordinates") + ylab("PredictProtein") + ggtitle("TBR1 Positive Controls")
p2 <- ggplot(tbr1.mr, aes(x = aapos.x, y = CADD.phred)) + geom_point(aes(colour=PC)) + xlab("Amino Acid Coordinates") + ylab("CADD Phred") + ggtitle("TBR1 Positive Controls")

# CADD vs PredictProtein
p3 <- ggplot(tbr1.mr, aes(x = PredictProtein, y = CADD.phred)) + geom_point(aes(colour=PC)) + xlab("PredictProtein") + ylab("CADD Phred") + ggtitle("TBR1 Positive Controls")
cor(tbr1.mr$CADD.phred, tbr1.mr$PredictProtein, method='spearman')

###
### Filter Negative controls
tbr1.mr %>% filter(ExonicFunc.refGene == "nonsynonymous SNV") -> tbr1.mr.nc
tbr1.mr.nc %>% mutate(NC = ifelse(
  (PredictProtein < 0 & CADD.phred < 15 & (exac03 > 0 | X1000g2012apr_all > 0 | esp6500si_all > 0) & LJB_MutationTaster < 0.99) | #Exac - occurring negative controls
    ((PredictProtein < -50 & CADD.phred < 5) & (tbr1.mr.nc$aapos.x > 200 & tbr1.mr.nc$aapos.x < 400) & LJB_MutationTaster < 0.99), # Novel negative controls within T-Box region
  yes="NegativeControl", no = FALSE)) -> tbr1.mr.nc
# tbr1.mr.nc %>% mutate(NC = ifelse(PredictProtein < 0 & CADD.phred < 10 ,yes="NegativeControl", no = FALSE)) -> tbr1.mr.nc

tbr1.mr.nc$NC[is.na(tbr1.mr.nc$NC)] <- FALSE

# Graph the Negative Controls
p4 <- ggplot(tbr1.mr.nc, aes(x = aapos.x, y = PredictProtein)) + geom_point(aes(colour=NC)) + xlab("Amino Acid Coordinates") + ylab("PredictProtein") + ggtitle("TBR1 Negative Controls")
p5 <- ggplot(tbr1.mr.nc, aes(x = aapos.x, y = CADD.phred)) + geom_point(aes(colour=NC)) + xlab("Amino Acid Coordinates") + ylab("CADD Phred") + ggtitle("TBR1 Negative Controls")
p6 <- ggplot(tbr1.mr.nc, aes(x = PredictProtein, y = CADD.phred)) + geom_point(aes(colour=NC)) + xlab("PredictProtein") + ylab("CADD Phred") + ggtitle("TBR1 Negative Controls")
###
###
###
# p8 <- ggplot(tbr1.mr, aes( x= as.numeric(gsub("TBR1:NM_[0-9]+:exon[0-9]+:c.[A-Z]([0-9]+)[A-Z]:p.[A-Z][0-9]+[A-Z]","\\1",tbr1.mr$AAChange.refGene,perl=TRUE)),
#                    y = as.numeric(CADD.phred))) + geom_point(aes(colour = PC)) + xlab("cDNA coords") + ylab("CADD Phred") + ggtitle("TBR1")

# ggplot(tbr1.mr, aes( x= as.numeric(gsub("TBR1:NM_[0-9]+:exon[0-9]+:c.[A-Z]([0-9]+)[A-Z]:p.[A-Z][0-9]+[A-Z]","\\1",tbr1.mr$AAChange.refGene,perl=TRUE)),
#                   y = as.numeric(PredictProtein))) + geom_point(aes(colour = PC)) + xlab("cDNA coords") + ylab("PredictProtein Score")+ ggtitle("TBR1")


p9 <- ggplot(tbr1.mr, aes( x= as.numeric(Start),y = as.numeric(CADD.phred))) + geom_point() + xlab('Genomic Coordinates') + ylab("CADD score") + 
  ggtitle("CADD scores by genomic coordinates")

p10 <- ggplot(tbr1.mr, aes(CADD.phred, PredictProtein)) + geom_point()

#What's that weird gap?
# tbr1.mr %>% select(c(1:4,11,17,18)) %>% slice(c(1500:2000))
cor(tbr1.mr$CADD.phred, tbr1.mr$PredictProtein, method='spearman')

tbr1.mr %>% filter(PredictProtein > 50) %>% filter(CADD.phred > 25) %>% select(-9)

# List coordinates of exons
v <- tbr1.mr$Start
split(unique(v),cumsum(c(1,diff(unique(v))!=1)))



# Output of tables
tbr1.mr.nc %>% filter(PC == "PositiveControl")  -> pcs
tbr1.mr.nc %>% filter(NC == "NegativeControl")  -> ncs
cons <- rbind(pcs,ncs)
write.table(cons, "tbr1controls.tsv", quote = FALSE, row.names=FALSE, sep = "\t")

tbr1.mr.nc %>% filter(aachange == "V356M" |aachange == "N374H" |aachange == "K228E" |aachange == "Q178E" | aachange == "Q418R" | aachange == "P542R") -> der
tbr1.mr.nc %>% filter(aachange == "V356M" |aachange == "N374H" |aachange == "K228E" |aachange == "Q178E" | aachange == "Q418R" | aachange == "P542R") %>% select(-c(11,12,13,14,16,17:26,28,29,31:33)) -> der
tbr1.mr.nc %>% filter(PC == "PositiveControl") %>% select(-c(8,11,12,13,14,16,17:26,28,29,31:33)) -> pcs
tbr1.mr.nc %>% filter(NC == "NegativeControl") %>% select(-c(8,11,12,13,14,16,17:26,28,29,31:33)) -> ncs
cons <- rbind(pcs,ncs)
write.table(cons, "tbr1controls_lite.tsv", quote = FALSE, row.names=FALSE, sep = "\t")
write.table(der, "deriziotis_vars_all_columns.tsv", quote = FALSE, row.names = FALSE, sep = "\t")
# pdf output
pdf(file="tbr1controls.pdf", height = 12, width = 17)
grid.arrange(tableGrob(pcs), top = "TBR1 Positive Controls" );p1;p2;p3;grid.arrange(tableGrob(ncs), top = "TBR1 Negative Controls");p4;p5;p6;p9
dev.off()


cor(tbr1.mr$LJB_MutationTaster, tbr1.mr$CADD.phred)


### IGV - get screenshots per variant -  is the alignment suspicious?
# Controls are the same region but in other subjects (use either controls or case subjects)
# Is there weird stuff going on in that region in controls???
# By next week
# Add parent DNA availabilty to aspire
# For example if some ms variants do not have parental dna might just skim them off, lower the cost, find other variants to fill out the list
# Ying - a lot of the LOF variants are splice sites - is this concerning?
# Disregarding dn status - do we have a lot more LOF per individual than other studies)
# Get information on CNV's - are they inherited etc?
# Look at UCSC - these splice sites - which exon are they in, what's the expression profile etc also expression of particular mRNAs

# Waiting on GSC information