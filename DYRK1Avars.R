##### Find variants for model organisms for DYRK1A

##### Setup
library(dplyr)
library(xtable)
library(gridExtra)
library(ggplot2)
setwd("/home/bcallaghan/NateDBCopy")
####

#### Load data
DYRK1A <- read.csv("/home/bcallaghan/DYRK1A_anno.hg19_multianno.csv")
DYRK1APP <- read.table("/home/bcallaghan/NateDBCopy/DYRK1APP.tsv")

natevars <- read.csv(file="outputs/DYRK1A/ndb_vars_DYRK1A.csv")
head(natevars)
natevars %>% select(c(1,2,3,4,5,6,7,9,10,11,36,41)) -> p 
png(filename="outputs/DYRK1A/natvars.png",width=1500, height = 26 * nrow(p))
grid.arrange(tableGrob(p))
dev.off()
####

#### Gene / Transcript Info (not used yet)
# Should add .bed file, alternative transcripts?
GENE == "DYRK1A" # 
TRANSCRIPT == "NM_000314" #Canonical Transcript
FASTA == NULL
BED == NULL
####

#### Fix annotations
x <- DYRK1A$Otherinfo
y <- strsplit(as.character(x), split = '\t')
y <- unlist(y)
y <- as.numeric(y[seq(2,length(y),2)])
DYRK1A$CADD.phred <- as.numeric(as.character(y))
####

##### Synonymous Mutations for Nate
natesyn <- c(167,92,131,101,107)
DYRK1A %>% 
  filter(ExonicFunc.refGene == "synonymous SNV") %>%
  filter(grepl("DYRK1A:NM_000314.+",AAChange.refGene)) -> DYRK1A.syn
DYRK1A.syn$aapos <- as.numeric(gsub(".*(DYRK1A:NM_000314:exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.[A-Z]([0-9]+)[A-Z],).+","\\2",
                              DYRK1A.syn$AAChange.refGene))
DYRK1A.syn %>% filter(aapos %in% natesyn) -> DYRK1A.syn
DYRK1A.syn$aachange <- gsub("DYRK1A:NM_000314:exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.([A-Z][0-9]+[A-Z]).+","\\1",DYRK1A.syn$AAChange.refGene,perl=TRUE)
#####

#### Filtering
DYRK1A$Func.refGene <- gsub("exonic;splicing","exonic",DYRK1A$Func.refGene)
# Multiple isoforms in DYRK1A so filter for canonical NM_000314
DYRK1A %>% filter(grepl("DYRK1A:NM_000314.+",AAChange.refGene)) -> DYRK1A.f
DYRK1A.f$AAChange.refGene
DYRK1A.f %>% mutate(aapos = as.numeric(gsub(".*(DYRK1A:NM_000314:exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.[A-Z]([0-9]+)[A-Z],).+","\\2",DYRK1A.f$AAChange.refGene))) -> DYRK1A.f
ggplot(DYRK1A, aes(x= aapos, y = CADD.phred)) + geom_point()
DYRK1A$aachange <- gsub("DYRK1A:NM_000314:exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.([A-Z][0-9]+[A-Z]).+","\\1",DYRK1A$AAChange.refGene,perl=TRUE)
# What's PredictProtein score look like in the context of a gene?
DYRK1APP %>% mutate(aapos = as.numeric(gsub('[A-Z]([0-9]+)[A-Z]','\\1',V1))) -> DYRK1APP
ggplot(DYRK1APP, aes(x = aapos, y=V2)) + geom_point()
####

#### Merge protein Predict and fix up new df
DYRK1APP$V1 <- as.character(DYRK1APP$V1)
DYRK1A$aachange
DYRK1A.mr <- merge(DYRK1A,DYRK1APP, by.x='aachange', by.y = 'V1', all=FALSE)
DYRK1A.mr$CADD.phred <- as.numeric(DYRK1A.mr$CADD.phred)
DYRK1A.mr$PredictProtein <- as.numeric(DYRK1A.mr$V2)
DYRK1A.mr$exac03 <- as.numeric(as.character(DYRK1A.mr$exac03))
DYRK1A.mr$exac03[is.na(DYRK1A.mr$exac03)] <- 0
DYRK1A.mr %>% arrange(Start) -> DYRK1A.mr
####

### Filtering 

###

##### Filter Positive Controls
DYRK1A.mr %>% mutate(PC = ifelse(PredictProtein > 75 & CADD.phred > 30 & (exac03 == 0 | is.na(exac03)) ,yes="PositiveControl", no = FALSE)) -> DYRK1A.mr
# Graph Predict Protein
ggplot(DYRK1A.mr, aes(x = aapos, y = PredictProtein)) + geom_point(aes()) + xlab("Amino Acid Coordinates") + ylab("PredictProtein") + ggtitle("DYRK1A")
# Graph CADD
ggplot(DYRK1A.mr, aes(x = aapos, y = CADD.phred)) + geom_point(aes()) + xlab("Amino Acid Coordinates") + ylab("CADD") + ggtitle("DYRK1A")
####

#### Graph the Positive Controls
p1 <- ggplot(DYRK1A.mr, aes(x = aapos, y = PredictProtein, colour = PC,shape = PC)) + 
  geom_point(size = 2.5) + xlab("Amino Acid Coordinates") + ylab("PredictProtein") + ggtitle("DYRK1A Positive Controls")
p2 <- ggplot(DYRK1A.mr, aes(x = aapos, y = CADD.phred, colour = PC,shape = PC)) + 
  geom_point(size = 2.5) + xlab("Amino Acid Coordinates") + ylab("CADD Phred") + ggtitle("DYRK1A Positive Controls")
# CADD vs PredictProtein
p3 <- ggplot(DYRK1A.mr, aes(x = PredictProtein, y = CADD.phred,colour = PC,shape = PC)) + 
  geom_point(size = 2.5) + xlab("PredictProtein") + ylab("CADD Phred") + ggtitle("DYRK1A Positive Controls")
cor(DYRK1A.mr$CADD.phred, DYRK1A.mr$PredictProtein, method='spearman')
####

#### Negative Controls Filtering
DYRK1A.mr %>% filter(ExonicFunc.refGene == "nonsynonymous SNV") -> DYRK1A.mr.nc
DYRK1A.mr.nc %>% mutate(NC = ifelse( test = ((PredictProtein < -50 & CADD.phred < 15 & exac03 > 0)|(PredictProtein < -50 & CADD.phred < 5)), 
                                   yes = "NegativeControl", no = FALSE)) -> DYRK1A.mr.nc
DYRK1A.mr.nc$NC[is.na(DYRK1A.mr.nc$NC)] <- FALSE
####

#### Graphing the Negative Controls
p4 <- ggplot(DYRK1A.mr.nc, aes(x = aapos, y = PredictProtein,colour=NC,shape = NC)) + 
  geom_point(size = 2.5) + xlab("Amino Acid Coordinates") + ylab("PredictProtein") + ggtitle("DYRK1A Negative Controls");p4
p5 <- ggplot(DYRK1A.mr.nc, aes(x = aapos, y = CADD.phred, colour = NC, shape = NC)) + 
  geom_point(size = 2.5) + xlab("Amino Acid Coordinates") + ylab("CADD Phred") + ggtitle("DYRK1A Negative Controls")
p6 <- ggplot(DYRK1A.mr.nc, aes(x = PredictProtein, y = CADD.phred, colour = NC, shape = NC)) + 
  geom_point(size = 2.5) + xlab("PredictProtein") + ylab("CADD Phred") + ggtitle("DYRK1A Negative Controls")
p9 <- ggplot(DYRK1A.mr, aes( x= as.numeric(Start),y = as.numeric(CADD.phred))) + 
  geom_point(size = 2.5) + xlab('Genomic Coordinates') + ylab("CADD score") + ggtitle("CADD by Genomic Coordinates")
p10 <- ggplot(DYRK1A.mr, aes(CADD.phred, PredictProtein)) + geom_point()

graph_variants <- function(anno.df,x.name,y.name,group,title){
  ggplot(anno.df, aes_string(x = x.name, y = y.name, colour = group, shape = group)) + 
    geom_point(size = 3) + 
    xlab(x.name) + ylab(y.name) + ggtitle(title)
}
graph_variants(DYRK1A,'aapos','CADD.phred','NC',"DYRK1A Negative Controls")

graph_variant_list <- function(anno.df,x.name,y.name,var.list,group,title, legend.title){
  anno.df %>% filter(Func.refGene == "exonic") -> anno.df.f
  anno.df.f %>% mutate(in.group = ifelse(aachange %in% var.list, yes = TRUE, no = FALSE)) -> anno.df.fm
  anno.df.fm$aapos <- as.numeric(gsub(".*(DYRK1A:NM_000314:exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.[A-Z]([0-9]+)[A-Z],).+","\\2",anno.df.fm$AAChange.refGene))
  anno.df.fm %>% filter(in.group == TRUE) -> anno.df.trus
  anno.df.trus %>% print()
  
#   group <- anno.df.fm$in.group
  ggplot(anno.df.fm, aes_string(x = x.name, y = y.name, colour = group)) + 
    geom_point(size = 3, alpha = 0.8) + 
  geom_point(aes_string(x = x.name, y = y.name, colour = group),data = anno.df.trus, size = 5) +
    xlab(x.name) + ylab(y.name) + ggtitle(title) + guides(title = "in.list") + 
#   theme(legend.title = element_text("asd"))+
  labs(colour = legend.title ) + guides(shape = FALSE) #+
  #scale_y_continuous( limits = c(0,50), expand = c(0,0) ) 
}
p11 <- graph_variant_list(DYRK1A.mr,'aapos','CADD.phred',DYRK1A.syn$aachange,"in.group","DYRK1A Synonymous Mutations","Synonymous");p11
p12 <- graph_variant_list(DYRK1A.mr,'aapos','PredictProtein',DYRK1A.syn$aachange,"in.group","DYRK1A Synonymous Mutations","Synonymous");p12

p11 <- graph_variant_list(DYRK1A.mr,'aapos','PredictProtein',pcs$aachange,"in.group","DYRK1A High Impact Mutations","High Impact");p11


DYRK1A.mr.nc %>% filter(PC == "PositiveControl") %>% select(aachange) -> pos.list
DYRK1A.mr.nc %>% filter(NC == "NegativeControl") %>% select(aachange) -> neg.list
DYRK1A %>% filter(ExonicFunc.refGene == "synonymous SNV") %>% select(aachange) -> syn.list

var.list <- c("L265L")
neg.list$aachange

# p11 <- graph_variant_list(DYRK1A,'aapos','CADD.phred',DYRK1A.syn$aachange,"in.group","DYRK1A Synonymous Mutations")
####Lit variants
### Filtering Nate Variants
natevars %>% mutate(stalt = paste0(Start,Ref,Alt)) -> natevars
DYRK1A.mr.nc %>% mutate(stalt = paste0(Start,Ref,Alt)) -> DYRK1A.mr.nc
DYRK1A.mr.nc %>% filter(stalt %in% natevars$stalt) -> DYRK1A.nate



#### Correlation / List coordinates of exons
cor(DYRK1A.mr$CADD.phred, DYRK1A.mr$PredictProtein, method='spearman')
v <- DYRK1A.mr$Start
split(unique(v),cumsum(c(1,diff(unique(v))!=1)))
####

#### Table output
DYRK1A.mr.nc %>% filter(PC == "PositiveControl")  -> pcs
DYRK1A.mr.nc %>% filter(NC == "NegativeControl")  -> ncs
cons <- rbind(pcs,ncs)
write.table(cons, "DYRK1Acontrols.tsv", quote = FALSE, row.names=FALSE, sep = "\t")
####

#### Controls Filtering
DYRK1A.mr.nc %>% filter(PC == "PositiveControl") %>% 
  select(aachange,Chr,Start,Ref,Alt,Func.refGene,ExonicFunc.refGene,snp135,exac03,CADD.phred,PredictProtein,PC,NC) -> pcs
DYRK1A.mr.nc %>% filter(NC == "NegativeControl") %>% 
  select(aachange,Chr,Start,Ref,Alt,Func.refGene,ExonicFunc.refGene,snp135,exac03,CADD.phred,PredictProtein,PC,NC) -> ncs

DYRK1A.syn$NC <- "SynonymousControl"
DYRK1A.syn$PC <- FALSE
# DYRK1A.syn$PredictProtein <- 0
DYRK1A.syn$PredictProtein <- DYRK1APP$V2[which(DYRK1APP$V1 %in% DYRK1A.syn$aachange)]
DYRK1A.syn

DYRK1A.nate %>% 
  select(aachange,Chr,Start,Ref,Alt,Func.refGene,ExonicFunc.refGene,snp135,exac03,CADD.phred,PredictProtein,PC,NC) -> lit

DYRK1A.syn %>% 
  select(aachange,Chr,Start,Ref,Alt,Func.refGene,ExonicFunc.refGene,snp135,exac03,CADD.phred,PredictProtein,PC,NC) -> scs
cons <- rbind(pcs,ncs,scs)

#### Table output
write.table(cons, "DYRK1Acontrols_lite.tsv", quote = FALSE, row.names=FALSE, sep = "\t")
write.table(lit,"DYRK1A_MARV_vars.tsv", quote = FALSE, row.names=FALSE, sep = "\t")
# write.table(der, "deriziotis_vars_all_columns.tsv", quote = FALSE, row.names = FALSE, sep = "\t")
####

#### Call graphing functions
# lit
p01 <- graph_variant_list(DYRK1A.mr,'aapos','CADD.phred',DYRK1A.nate$aachange,"in.group","DYRK1A Literature Mutations","Literature");
p02 <- graph_variant_list(DYRK1A.mr,'aapos','PredictProtein',DYRK1A.nate$aachange,"in.group","DYRK1A Literature Mutations","Literature");
p03 <- graph_variant_list(DYRK1A.mr,'PredictProtein','CADD.phred',DYRK1A.nate$aachange,"in.group","DYRK1A Literature Mutations","Literature");
#syn
p1 <- graph_variant_list(DYRK1A.mr,'aapos','CADD.phred',DYRK1A.syn$aachange,"in.group","DYRK1A Synonymous Mutations","Synonymous");
p2 <- graph_variant_list(DYRK1A.mr,'aapos','PredictProtein',DYRK1A.syn$aachange,"in.group","DYRK1A Synonymous Mutations","Synonymous");
p3 <- graph_variant_list(DYRK1A.mr,'PredictProtein','CADD.phred',DYRK1A.syn$aachange,"in.group","DYRK1A Synonymous Mutations","Synonymous");

#pcs
p4 <- graph_variant_list(DYRK1A.mr,'aapos','CADD.phred',pcs$aachange,"in.group","DYRK1A High Impact Mutations","High Impact");
p5 <- graph_variant_list(DYRK1A.mr,'aapos','PredictProtein',pcs$aachange,"in.group","DYRK1A High Impact Mutations","High Impact");
p6 <- graph_variant_list(DYRK1A.mr,'PredictProtein','CADD.phred',pcs$aachange,"in.group","DYRK1A High Impact Mutations","High Impact");

#ncs
p7 <- graph_variant_list(DYRK1A.mr,'aapos','CADD.phred',ncs$aachange,"in.group","DYRK1A Low Impact Mutations","Low Impact");
p8 <- graph_variant_list(DYRK1A.mr,'aapos','PredictProtein',ncs$aachange,"in.group","DYRK1A Low Impact Mutations","Low Impact");
p9 <- graph_variant_list(DYRK1A.mr,'PredictProtein','CADD.phred',ncs$aachange,"in.group","Low Impact Mutations","Low Impact");
####

## PDF output 2
pdf(file="outputs/DYRK1A/DYRK1Acontrols3.pdf", height = 12, width = 17)
grid.arrange(tableGrob(lit), top = "DYRK1A Literature Mutations" );p01;p02;p03
grid.arrange(tableGrob(scs), top = "DYRK1A Synonymous Mutations" );p1;p2;p3
grid.arrange(tableGrob(pcs), top = "DYRK1A High Impact Mutations" );p4;p5;p6
grid.arrange(tableGrob(ncs), top = "DYRK1A Low Impact Mutations");p7;p8;p9;p9
#p10
dev.off()


#### PDF output
# pdf(file="outputs/DYRK1A/DYRK1Acontrols.pdf", height = 12, width = 17)
# grid.arrange(tableGrob(scs), top = "DYRK1A Synonymous Controls" );p11;p12
# grid.arrange(tableGrob(pcs), top = "DYRK1A Positive Controls" );p1;p2;p3
# grid.arrange(tableGrob(ncs), top = "DYRK1A Negative Controls");p4;p5;p6;p9
# #p10
# dev.off()
####

#### What happens when we compare same aa changes (different nuc change)?
DYRK1A.mr.nc %>% filter(duplicated(aachange,fromLast=FALSE) | duplicated(aachange,fromLast=TRUE)) %>%
  arrange(aapos.x)-> sameaa
ggplot(sameaa, aes(x = aapos.x, y = CADD.phred)) + geom_point(aes(colour=PredictProtein)) + 
  xlab("AA Pos") + ylab("CADD Phred") + ggtitle("DYRK1A Same aa vars") + 
  scale_color_gradient2(low = "green", mid = "yellow", high = "red" )
####

# cor(DYRK1A.mr$LJB_MutationTaster, DYRK1A.mr$CADD.phred)

#### Notes
# IGV - get screenshots per variant -  is the alignment suspicious?
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




