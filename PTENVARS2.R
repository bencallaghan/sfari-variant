##### Find variants for model organisms for PTEN

##### Setup
library(dplyr)
library(xtable)
library(gridExtra)
library(ggplot2)
setwd("/home/bcallaghan/NateDBCopy")
####

#### Load data
PTEN <- read.csv("/home/bcallaghan/pten_anno.hg19_multianno.csv")
PTENPP <- read.table("/home/bcallaghan/NateDBCopy/PTENPP.tsv")

natevars <- read.csv(file="outputs/PTEN/ndb_vars_PTEN.csv")
head(natevars)
natevars %>% select(c(1,2,3,4,5,6,7,9,10,11,36,41)) -> p 
png(filename="outputs/PTEN/natvars.png",width=1500, height = 26 * nrow(p))
grid.arrange(tableGrob(p))
dev.off()
####

#### Gene / Transcript Info (not used yet)
# Should add .bed file, alternative transcripts?
GENE == "PTEN" # 
TRANSCRIPT == "NM_000314" #Canonical Transcript
FASTA == NULL
BED == NULL
####

#### Fix annotations
x <- PTEN$Otherinfo
y <- strsplit(as.character(x), split = '\t')
y <- unlist(y)
y <- as.numeric(y[seq(2,length(y),2)])
PTEN$CADD.phred <- as.numeric(as.character(y))
####

##### Synonymous Mutations for Nate
natesyn <- c(167,92,131,101,107)
PTEN %>% 
  filter(ExonicFunc.refGene == "synonymous SNV") %>%
  filter(grepl("PTEN:NM_000314.+",AAChange.refGene)) -> PTEN.syn
PTEN.syn$aapos <- as.numeric(gsub(".*(PTEN:NM_000314:exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.[A-Z]([0-9]+)[A-Z],).+","\\2",
                              PTEN.syn$AAChange.refGene))
PTEN.syn %>% filter(aapos %in% natesyn) -> PTEN.syn
PTEN.syn$aachange <- gsub("PTEN:NM_000314:exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.([A-Z][0-9]+[A-Z]).+","\\1",PTEN.syn$AAChange.refGene,perl=TRUE)
#####

#### Filtering
PTEN$Func.refGene <- gsub("exonic;splicing","exonic",PTEN$Func.refGene)
# Multiple isoforms in PTEN so filter for canonical NM_000314
PTEN %>% filter(grepl("PTEN:NM_000314.+",AAChange.refGene)) -> PTEN.f
PTEN.f$AAChange.refGene
PTEN.f %>% mutate(aapos = as.numeric(gsub(".*(PTEN:NM_000314:exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.[A-Z]([0-9]+)[A-Z],).+","\\2",PTEN.f$AAChange.refGene))) -> PTEN.f
ggplot(PTEN, aes(x= aapos, y = CADD.phred)) + geom_point()
PTEN$aachange <- gsub("PTEN:NM_000314:exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.([A-Z][0-9]+[A-Z]).+","\\1",PTEN$AAChange.refGene,perl=TRUE)
# What's PredictProtein score look like in the context of a gene?
PTENPP %>% mutate(aapos = as.numeric(gsub('[A-Z]([0-9]+)[A-Z]','\\1',V1))) -> PTENPP
ggplot(PTENPP, aes(x = aapos, y=V2)) + geom_point()
####

#### Merge protein Predict and fix up new df
PTENPP$V1 <- as.character(PTENPP$V1)
PTEN$aachange
PTEN.mr <- merge(PTEN,PTENPP, by.x='aachange', by.y = 'V1', all=FALSE)
PTEN.mr$CADD.phred <- as.numeric(PTEN.mr$CADD.phred)
PTEN.mr$PredictProtein <- as.numeric(PTEN.mr$V2)
PTEN.mr$exac03 <- as.numeric(as.character(PTEN.mr$exac03))
PTEN.mr$exac03[is.na(PTEN.mr$exac03)] <- 0
PTEN.mr %>% arrange(Start) -> PTEN.mr
####

### Filtering 

###

##### Filter Positive Controls
PTEN.mr %>% mutate(PC = ifelse(PredictProtein > 75 & CADD.phred > 30 & (exac03 == 0 | is.na(exac03)) ,yes="PositiveControl", no = FALSE)) -> PTEN.mr
# Graph Predict Protein
ggplot(PTEN.mr, aes(x = aapos, y = PredictProtein)) + geom_point(aes()) + xlab("Amino Acid Coordinates") + ylab("PredictProtein") + ggtitle("PTEN")
# Graph CADD
ggplot(PTEN.mr, aes(x = aapos, y = CADD.phred)) + geom_point(aes()) + xlab("Amino Acid Coordinates") + ylab("CADD") + ggtitle("PTEN")
####

#### Graph the Positive Controls
p1 <- ggplot(PTEN.mr, aes(x = aapos, y = PredictProtein, colour = PC,shape = PC)) + 
  geom_point(size = 2.5) + xlab("Amino Acid Coordinates") + ylab("PredictProtein") + ggtitle("PTEN Positive Controls")
p2 <- ggplot(PTEN.mr, aes(x = aapos, y = CADD.phred, colour = PC,shape = PC)) + 
  geom_point(size = 2.5) + xlab("Amino Acid Coordinates") + ylab("CADD Phred") + ggtitle("PTEN Positive Controls")
# CADD vs PredictProtein
p3 <- ggplot(PTEN.mr, aes(x = PredictProtein, y = CADD.phred,colour = PC,shape = PC)) + 
  geom_point(size = 2.5) + xlab("PredictProtein") + ylab("CADD Phred") + ggtitle("PTEN Positive Controls")
cor(PTEN.mr$CADD.phred, PTEN.mr$PredictProtein, method='spearman')
####

#### Negative Controls Filtering
PTEN.mr %>% filter(ExonicFunc.refGene == "nonsynonymous SNV") -> PTEN.mr.nc
PTEN.mr.nc %>% mutate(NC = ifelse( test = ((PredictProtein < -50 & CADD.phred < 15 & exac03 > 0)|(PredictProtein < -50 & CADD.phred < 5)), 
                                   yes = "NegativeControl", no = FALSE)) -> PTEN.mr.nc
PTEN.mr.nc$NC[is.na(PTEN.mr.nc$NC)] <- FALSE
####

#### Graphing the Negative Controls
p4 <- ggplot(PTEN.mr.nc, aes(x = aapos, y = PredictProtein,colour=NC,shape = NC)) + 
  geom_point(size = 2.5) + xlab("Amino Acid Coordinates") + ylab("PredictProtein") + ggtitle("PTEN Negative Controls");p4
p5 <- ggplot(PTEN.mr.nc, aes(x = aapos, y = CADD.phred, colour = NC, shape = NC)) + 
  geom_point(size = 2.5) + xlab("Amino Acid Coordinates") + ylab("CADD Phred") + ggtitle("PTEN Negative Controls")
p6 <- ggplot(PTEN.mr.nc, aes(x = PredictProtein, y = CADD.phred, colour = NC, shape = NC)) + 
  geom_point(size = 2.5) + xlab("PredictProtein") + ylab("CADD Phred") + ggtitle("PTEN Negative Controls")
p9 <- ggplot(PTEN.mr, aes( x= as.numeric(Start),y = as.numeric(CADD.phred))) + 
  geom_point(size = 2.5) + xlab('Genomic Coordinates') + ylab("CADD score") + ggtitle("CADD by Genomic Coordinates")
p10 <- ggplot(PTEN.mr, aes(CADD.phred, PredictProtein)) + geom_point()

graph_variants <- function(anno.df,x.name,y.name,group,title){
  ggplot(anno.df, aes_string(x = x.name, y = y.name, colour = group, shape = group)) + 
    geom_point(size = 3) + 
    xlab(x.name) + ylab(y.name) + ggtitle(title)
}
graph_variants(PTEN,'aapos','CADD.phred','NC',"PTEN Negative Controls")

graph_variant_list <- function(anno.df,x.name,y.name,var.list,group,title, legend.title){
  anno.df %>% filter(Func.refGene == "exonic") -> anno.df.f
  anno.df.f %>% mutate(in.group = ifelse(aachange %in% var.list, yes = TRUE, no = FALSE)) -> anno.df.fm
  anno.df.fm$aapos <- as.numeric(gsub(".*(PTEN:NM_000314:exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.[A-Z]([0-9]+)[A-Z],).+","\\2",anno.df.fm$AAChange.refGene))
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
p11 <- graph_variant_list(PTEN.mr,'aapos','CADD.phred',PTEN.syn$aachange,"in.group","PTEN Synonymous Mutations","Synonymous");p11
p12 <- graph_variant_list(PTEN.mr,'aapos','PredictProtein',PTEN.syn$aachange,"in.group","PTEN Synonymous Mutations","Synonymous");p12

p11 <- graph_variant_list(PTEN.mr,'aapos','PredictProtein',pcs$aachange,"in.group","PTEN High Impact Mutations","High Impact");p11


PTEN.mr.nc %>% filter(PC == "PositiveControl") %>% select(aachange) -> pos.list
PTEN.mr.nc %>% filter(NC == "NegativeControl") %>% select(aachange) -> neg.list
PTEN %>% filter(ExonicFunc.refGene == "synonymous SNV") %>% select(aachange) -> syn.list

var.list <- c("L265L")
neg.list$aachange

# p11 <- graph_variant_list(PTEN,'aapos','CADD.phred',PTEN.syn$aachange,"in.group","PTEN Synonymous Mutations")
####Lit variants
### Filtering Nate Variants
natevars %>% mutate(stalt = paste0(Start,Ref,Alt)) -> natevars
PTEN.mr.nc %>% mutate(stalt = paste0(Start,Ref,Alt)) -> PTEN.mr.nc
PTEN.mr.nc %>% filter(stalt %in% natevars$stalt) -> PTEN.nate



#### Correlation / List coordinates of exons
cor(PTEN.mr$CADD.phred, PTEN.mr$PredictProtein, method='spearman')
v <- PTEN.mr$Start
split(unique(v),cumsum(c(1,diff(unique(v))!=1)))
####

#### Table output
PTEN.mr.nc %>% filter(PC == "PositiveControl")  -> pcs
PTEN.mr.nc %>% filter(NC == "NegativeControl")  -> ncs
cons <- rbind(pcs,ncs)
write.table(cons, "PTENcontrols.tsv", quote = FALSE, row.names=FALSE, sep = "\t")
####

#### Controls Filtering
PTEN.mr.nc %>% filter(PC == "PositiveControl") %>% 
  select(aachange,Chr,Start,Ref,Alt,Func.refGene,ExonicFunc.refGene,snp135,exac03,CADD.phred,PredictProtein,PC,NC) -> pcs
PTEN.mr.nc %>% filter(NC == "NegativeControl") %>% 
  select(aachange,Chr,Start,Ref,Alt,Func.refGene,ExonicFunc.refGene,snp135,exac03,CADD.phred,PredictProtein,PC,NC) -> ncs

PTEN.syn$NC <- "SynonymousControl"
PTEN.syn$PC <- FALSE
# PTEN.syn$PredictProtein <- 0
PTEN.syn$PredictProtein <- PTENPP$V2[which(PTENPP$V1 %in% PTEN.syn$aachange)]
PTEN.syn

PTEN.nate %>% 
  select(aachange,Chr,Start,Ref,Alt,Func.refGene,ExonicFunc.refGene,snp135,exac03,CADD.phred,PredictProtein,PC,NC) -> lit

PTEN.syn %>% 
  select(aachange,Chr,Start,Ref,Alt,Func.refGene,ExonicFunc.refGene,snp135,exac03,CADD.phred,PredictProtein,PC,NC) -> scs
cons <- rbind(pcs,ncs,scs)

#### Table output
write.table(cons, "PTENcontrols_lite.tsv", quote = FALSE, row.names=FALSE, sep = "\t")
write.table(lit,"PTEN_MARV_vars.tsv", quote = FALSE, row.names=FALSE, sep = "\t")
# write.table(der, "deriziotis_vars_all_columns.tsv", quote = FALSE, row.names = FALSE, sep = "\t")
####

#### Call graphing functions
# lit
p01 <- graph_variant_list(PTEN.mr,'aapos','CADD.phred',PTEN.nate$aachange,"in.group","PTEN Literature Mutations","Literature");
p02 <- graph_variant_list(PTEN.mr,'aapos','PredictProtein',PTEN.nate$aachange,"in.group","PTEN Literature Mutations","Literature");
p03 <- graph_variant_list(PTEN.mr,'PredictProtein','CADD.phred',PTEN.nate$aachange,"in.group","PTEN Literature Mutations","Literature");
#syn
p1 <- graph_variant_list(PTEN.mr,'aapos','CADD.phred',PTEN.syn$aachange,"in.group","PTEN Synonymous Mutations","Synonymous");
p2 <- graph_variant_list(PTEN.mr,'aapos','PredictProtein',PTEN.syn$aachange,"in.group","PTEN Synonymous Mutations","Synonymous");
p3 <- graph_variant_list(PTEN.mr,'PredictProtein','CADD.phred',PTEN.syn$aachange,"in.group","PTEN Synonymous Mutations","Synonymous");

#pcs
p4 <- graph_variant_list(PTEN.mr,'aapos','CADD.phred',pcs$aachange,"in.group","PTEN High Impact Mutations","High Impact");
p5 <- graph_variant_list(PTEN.mr,'aapos','PredictProtein',pcs$aachange,"in.group","PTEN High Impact Mutations","High Impact");
p6 <- graph_variant_list(PTEN.mr,'PredictProtein','CADD.phred',pcs$aachange,"in.group","PTEN High Impact Mutations","High Impact");

#ncs
p7 <- graph_variant_list(PTEN.mr,'aapos','CADD.phred',ncs$aachange,"in.group","PTEN Low Impact Mutations","Low Impact");
p8 <- graph_variant_list(PTEN.mr,'aapos','PredictProtein',ncs$aachange,"in.group","PTEN Low Impact Mutations","Low Impact");
p9 <- graph_variant_list(PTEN.mr,'PredictProtein','CADD.phred',ncs$aachange,"in.group","Low Impact Mutations","Low Impact");
####

## PDF output 2
pdf(file="outputs/PTEN/PTENcontrols3.pdf", height = 12, width = 17)
grid.arrange(tableGrob(lit), top = "PTEN Literature Mutations" );p01;p02;p03
grid.arrange(tableGrob(scs), top = "PTEN Synonymous Mutations" );p1;p2;p3
grid.arrange(tableGrob(pcs), top = "PTEN High Impact Mutations" );p4;p5;p6
grid.arrange(tableGrob(ncs), top = "PTEN Low Impact Mutations");p7;p8;p9;p9
#p10
dev.off()


#### PDF output
# pdf(file="outputs/PTEN/PTENcontrols.pdf", height = 12, width = 17)
# grid.arrange(tableGrob(scs), top = "PTEN Synonymous Controls" );p11;p12
# grid.arrange(tableGrob(pcs), top = "PTEN Positive Controls" );p1;p2;p3
# grid.arrange(tableGrob(ncs), top = "PTEN Negative Controls");p4;p5;p6;p9
# #p10
# dev.off()
####

#### What happens when we compare same aa changes (different nuc change)?
PTEN.mr.nc %>% filter(duplicated(aachange,fromLast=FALSE) | duplicated(aachange,fromLast=TRUE)) %>%
  arrange(aapos.x)-> sameaa
ggplot(sameaa, aes(x = aapos.x, y = CADD.phred)) + geom_point(aes(colour=PredictProtein)) + 
  xlab("AA Pos") + ylab("CADD Phred") + ggtitle("PTEN Same aa vars") + 
  scale_color_gradient2(low = "green", mid = "yellow", high = "red" )
####

# cor(PTEN.mr$LJB_MutationTaster, PTEN.mr$CADD.phred)

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




