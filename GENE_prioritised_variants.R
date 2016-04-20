##### Find variants for model organisms for GENE.vars

##### Setup
library(dplyr)
library(xtable)
library(gridExtra)
library(ggplot2)
setwd("/home/bcallaghan/NateDBCopy")
####

 

#### Load data
GENE <- "SYNGAP1"
CHROM <- 6
TRANSCRIPT <-  "NM_006772"
INPUTDIR <- paste0("./inputs/",GENE,"/")
getwd()
FASTA == NULL
BED == NULL

anno.path <- paste0(INPUTDIR,GENE,"_anno.hg19_multianno.csv")
GENE.vars <- read.csv(anno.path)
pp.path <- paste0(INPUTDIR,GENE,"PP.tsv")
GENE.PP <- read.table(pp.path)

nate.path <- paste0(INPUTDIR,GENE,".hg19_multianno.csv")
natevars <- read.csv(nate.path)
head(natevars)
natevars %>% select(c(1,2,3,4,5,6,7,9,10,11,36,41)) -> p 
png(filename="outputs/natvars.png",width=1500, height = 26 * nrow(p))
grid.arrange(tableGrob(p))
dev.off()
####

#### Gene / Transcript Info (not used yet)
# # Should add .bed file, alternative transcripts?
# GENE == "GENE.vars" # 
 #Canonical Transcript

####

#### Fix annotations
x <- GENE.vars$Otherinfo
y <- strsplit(as.character(x), split = '\t')
y <- unlist(y)
y <- as.numeric(y[seq(2,length(y),2)])
GENE.vars$CADD.phred <- as.numeric(as.character(y))
####

##### Synonymous Mutations for Nate
natesyn <- natevars$Start
GENE.vars %>% 
  filter(ExonicFunc.refGene == "synonymous SNV") %>%
  filter(grepl("GENE.vars:NM_000314.+",AAChange.refGene)) -> GENE.vars.syn
GENE.vars.syn$aapos <- as.numeric(gsub(".*(GENE.vars:NM_000314:exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.[A-Z]([0-9]+)[A-Z],).+","\\2",
                              GENE.vars.syn$AAChange.refGene))
GENE.vars.syn %>% filter(aapos %in% natesyn) -> GENE.vars.syn
GENE.vars.syn$aachange <- gsub("GENE.vars:NM_000314:exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.([A-Z][0-9]+[A-Z]).+","\\1",GENE.vars.syn$AAChange.refGene,perl=TRUE)
#####

#### Filtering
GENE.vars$Func.refGene <- gsub("exonic;splicing","exonic",GENE.vars$Func.refGene)
# Multiple isoforms in GENE.vars so filter for canonical NM_000314
GENE.vars %>% filter(grepl(paste0(GENE,":",TRANSCRIPT,".+"),AAChange.refGene)) -> GENE.vars.f
GENE.vars.f$AAChange.refGene
# GENE.vars.f %>% mutate(aapos = as.numeric(gsub(".*(GENE.vars:NM_000314:exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.[A-Z]([0-9]+)[A-Z],).+","\\2",GENE.vars.f$AAChange.refGene))) -> GENE.vars.f

#********************
GENE.vars.f$aapos <- as.numeric(gsub(".*GENE.vars:NM_000314:exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.[A-Z]([0-9]+)[A-Z],).+","\\2",GENE.vars.f$AAChange.refGene))
GENE.vars.f$aapos <- as.numeric(gsub(paste0( ".*" , GENE , ":" , TRANSCRIPT , ":exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.[A-Z]([0-9]+)[A-Z],).+"),"\\2",GENE.vars.f$AAChange.refGene)) 
# *********************

ggplot(GENE.vars.f, aes(x= aapos, y = CADD.phred)) + geom_point()
GENE.vars$aachange <- gsub("GENE.vars:NM_000314:exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.([A-Z][0-9]+[A-Z]).+","\\1",GENE.vars$AAChange.refGene,perl=TRUE)
# What's PredictProtein score look like in the context of a gene?
GENE.PP %>% mutate(aapos = as.numeric(gsub('[A-Z]([0-9]+)[A-Z]','\\1',V1))) -> GENE.PP
ggplot(GENE.PP, aes(x = aapos, y=V2)) + geom_point()
####

#### Merge protein Predict and fix up new df
GENE.PP$V1 <- as.character(GENE.PP$V1)
GENE.vars$aachange
GENE.vars.mr <- merge(GENE.vars,GENE.PP, by.x='aachange', by.y = 'V1', all=FALSE)
GENE.vars.mr$CADD.phred <- as.numeric(GENE.vars.mr$CADD.phred)
GENE.vars.mr$PredictProtein <- as.numeric(GENE.vars.mr$V2)
GENE.vars.mr$exac03 <- as.numeric(as.character(GENE.vars.mr$exac03))
GENE.vars.mr$exac03[is.na(GENE.vars.mr$exac03)] <- 0
GENE.vars.mr %>% arrange(Start) -> GENE.vars.mr
####

### Filtering 

###

##### Filter Positive Controls
GENE.vars.mr %>% mutate(PC = ifelse(PredictProtein > 75 & CADD.phred > 30 & (exac03 == 0 | is.na(exac03)) ,yes="PositiveControl", no = FALSE)) -> GENE.vars.mr
# Graph Predict Protein
ggplot(GENE.vars.mr, aes(x = aapos, y = PredictProtein)) + geom_point(aes()) + xlab("Amino Acid Coordinates") + ylab("PredictProtein") + ggtitle("GENE.vars")
# Graph CADD
ggplot(GENE.vars.mr, aes(x = aapos, y = CADD.phred)) + geom_point(aes()) + xlab("Amino Acid Coordinates") + ylab("CADD") + ggtitle("GENE.vars")
####

#### Graph the Positive Controls
p1 <- ggplot(GENE.vars.mr, aes(x = aapos, y = PredictProtein, colour = PC,shape = PC)) + 
  geom_point(size = 2.5) + xlab("Amino Acid Coordinates") + ylab("PredictProtein") + ggtitle("GENE.vars Positive Controls")
p2 <- ggplot(GENE.vars.mr, aes(x = aapos, y = CADD.phred, colour = PC,shape = PC)) + 
  geom_point(size = 2.5) + xlab("Amino Acid Coordinates") + ylab("CADD Phred") + ggtitle("GENE.vars Positive Controls")
# CADD vs PredictProtein
p3 <- ggplot(GENE.vars.mr, aes(x = PredictProtein, y = CADD.phred,colour = PC,shape = PC)) + 
  geom_point(size = 2.5) + xlab("PredictProtein") + ylab("CADD Phred") + ggtitle("GENE.vars Positive Controls")
cor(GENE.vars.mr$CADD.phred, GENE.vars.mr$PredictProtein, method='spearman')
####

#### Negative Controls Filtering
GENE.vars.mr %>% filter(ExonicFunc.refGene == "nonsynonymous SNV") -> GENE.vars.mr.nc
GENE.vars.mr.nc %>% mutate(NC = ifelse( test = ((PredictProtein < -50 & CADD.phred < 15 & exac03 > 0)|(PredictProtein < -50 & CADD.phred < 5)), 
                                   yes = "NegativeControl", no = FALSE)) -> GENE.vars.mr.nc
GENE.vars.mr.nc$NC[is.na(GENE.vars.mr.nc$NC)] <- FALSE
####

#### Graphing the Negative Controls
p4 <- ggplot(GENE.vars.mr.nc, aes(x = aapos, y = PredictProtein,colour=NC,shape = NC)) + 
  geom_point(size = 2.5) + xlab("Amino Acid Coordinates") + ylab("PredictProtein") + ggtitle("GENE.vars Negative Controls");p4
p5 <- ggplot(GENE.vars.mr.nc, aes(x = aapos, y = CADD.phred, colour = NC, shape = NC)) + 
  geom_point(size = 2.5) + xlab("Amino Acid Coordinates") + ylab("CADD Phred") + ggtitle("GENE.vars Negative Controls")
p6 <- ggplot(GENE.vars.mr.nc, aes(x = PredictProtein, y = CADD.phred, colour = NC, shape = NC)) + 
  geom_point(size = 2.5) + xlab("PredictProtein") + ylab("CADD Phred") + ggtitle("GENE.vars Negative Controls")
p9 <- ggplot(GENE.vars.mr, aes( x= as.numeric(Start),y = as.numeric(CADD.phred))) + 
  geom_point(size = 2.5) + xlab('Genomic Coordinates') + ylab("CADD score") + ggtitle("CADD by Genomic Coordinates")
p10 <- ggplot(GENE.vars.mr, aes(CADD.phred, PredictProtein)) + geom_point()

graph_variants <- function(anno.df,x.name,y.name,group,title){
  ggplot(anno.df, aes_string(x = x.name, y = y.name, colour = group, shape = group)) + 
    geom_point(size = 3) + 
    xlab(x.name) + ylab(y.name) + ggtitle(title)
}
graph_variants(GENE.vars,'aapos','CADD.phred','NC',"GENE.vars Negative Controls")

graph_variant_list <- function(anno.df,x.name,y.name,var.list,group,title, legend.title){
  anno.df %>% filter(Func.refGene == "exonic") -> anno.df.f
  anno.df.f %>% mutate(in.group = ifelse(aachange %in% var.list, yes = TRUE, no = FALSE)) -> anno.df.fm
  anno.df.fm$aapos <- as.numeric(gsub(".*(GENE.vars:NM_000314:exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.[A-Z]([0-9]+)[A-Z],).+","\\2",anno.df.fm$AAChange.refGene))
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
p11 <- graph_variant_list(GENE.vars.mr,'aapos','CADD.phred',GENE.vars.syn$aachange,"in.group","GENE.vars Synonymous Mutations","Synonymous");p11
p12 <- graph_variant_list(GENE.vars.mr,'aapos','PredictProtein',GENE.vars.syn$aachange,"in.group","GENE.vars Synonymous Mutations","Synonymous");p12

p11 <- graph_variant_list(GENE.vars.mr,'aapos','PredictProtein',pcs$aachange,"in.group","GENE.vars High Impact Mutations","High Impact");p11


GENE.vars.mr.nc %>% filter(PC == "PositiveControl") %>% select(aachange) -> pos.list
GENE.vars.mr.nc %>% filter(NC == "NegativeControl") %>% select(aachange) -> neg.list
GENE.vars %>% filter(ExonicFunc.refGene == "synonymous SNV") %>% select(aachange) -> syn.list

var.list <- c("L265L")
neg.list$aachange

# p11 <- graph_variant_list(GENE.vars,'aapos','CADD.phred',GENE.vars.syn$aachange,"in.group","GENE.vars Synonymous Mutations")
####Lit variants
### Filtering Nate Variants
natevars %>% mutate(stalt = paste0(Start,Ref,Alt)) -> natevars
GENE.vars.mr.nc %>% mutate(stalt = paste0(Start,Ref,Alt)) -> GENE.vars.mr.nc
GENE.vars.mr.nc %>% filter(stalt %in% natevars$stalt) -> GENE.vars.nate



#### Correlation / List coordinates of exons
cor(GENE.vars.mr$CADD.phred, GENE.vars.mr$PredictProtein, method='spearman')
v <- GENE.vars.mr$Start
split(unique(v),cumsum(c(1,diff(unique(v))!=1)))
####

#### Table output
GENE.vars.mr.nc %>% filter(PC == "PositiveControl")  -> pcs
GENE.vars.mr.nc %>% filter(NC == "NegativeControl")  -> ncs
cons <- rbind(pcs,ncs)
write.table(cons, "GENE.varscontrols.tsv", quote = FALSE, row.names=FALSE, sep = "\t")
####

#### Controls Filtering
GENE.vars.mr.nc %>% filter(PC == "PositiveControl") %>% 
  select(aachange,Chr,Start,Ref,Alt,Func.refGene,ExonicFunc.refGene,snp135,exac03,CADD.phred,PredictProtein,PC,NC) -> pcs
GENE.vars.mr.nc %>% filter(NC == "NegativeControl") %>% 
  select(aachange,Chr,Start,Ref,Alt,Func.refGene,ExonicFunc.refGene,snp135,exac03,CADD.phred,PredictProtein,PC,NC) -> ncs

GENE.vars.syn$NC <- "SynonymousControl"
GENE.vars.syn$PC <- FALSE
# GENE.vars.syn$PredictProtein <- 0
GENE.vars.syn$PredictProtein <- GENE.PP$V2[which(GENE.PP$V1 %in% GENE.vars.syn$aachange)]
GENE.vars.syn

GENE.vars.nate %>% 
  select(aachange,Chr,Start,Ref,Alt,Func.refGene,ExonicFunc.refGene,snp135,exac03,CADD.phred,PredictProtein,PC,NC) -> lit

GENE.vars.syn %>% 
  select(aachange,Chr,Start,Ref,Alt,Func.refGene,ExonicFunc.refGene,snp135,exac03,CADD.phred,PredictProtein,PC,NC) -> scs
cons <- rbind(pcs,ncs,scs)

#### Table output
write.table(cons, "GENE.varscontrols_lite.tsv", quote = FALSE, row.names=FALSE, sep = "\t")
write.table(lit,"GENE.vars_MARV_vars.tsv", quote = FALSE, row.names=FALSE, sep = "\t")
# write.table(der, "deriziotis_vars_all_columns.tsv", quote = FALSE, row.names = FALSE, sep = "\t")
####

#### Call graphing functions
# lit
p01 <- graph_variant_list(GENE.vars.mr,'aapos','CADD.phred',GENE.vars.nate$aachange,"in.group","GENE.vars Literature Mutations","Literature");
p02 <- graph_variant_list(GENE.vars.mr,'aapos','PredictProtein',GENE.vars.nate$aachange,"in.group","GENE.vars Literature Mutations","Literature");
p03 <- graph_variant_list(GENE.vars.mr,'PredictProtein','CADD.phred',GENE.vars.nate$aachange,"in.group","GENE.vars Literature Mutations","Literature");
#syn
p1 <- graph_variant_list(GENE.vars.mr,'aapos','CADD.phred',GENE.vars.syn$aachange,"in.group","GENE.vars Synonymous Mutations","Synonymous");
p2 <- graph_variant_list(GENE.vars.mr,'aapos','PredictProtein',GENE.vars.syn$aachange,"in.group","GENE.vars Synonymous Mutations","Synonymous");
p3 <- graph_variant_list(GENE.vars.mr,'PredictProtein','CADD.phred',GENE.vars.syn$aachange,"in.group","GENE.vars Synonymous Mutations","Synonymous");

#pcs
p4 <- graph_variant_list(GENE.vars.mr,'aapos','CADD.phred',pcs$aachange,"in.group","GENE.vars High Impact Mutations","High Impact");
p5 <- graph_variant_list(GENE.vars.mr,'aapos','PredictProtein',pcs$aachange,"in.group","GENE.vars High Impact Mutations","High Impact");
p6 <- graph_variant_list(GENE.vars.mr,'PredictProtein','CADD.phred',pcs$aachange,"in.group","GENE.vars High Impact Mutations","High Impact");

#ncs
p7 <- graph_variant_list(GENE.vars.mr,'aapos','CADD.phred',ncs$aachange,"in.group","GENE.vars Low Impact Mutations","Low Impact");
p8 <- graph_variant_list(GENE.vars.mr,'aapos','PredictProtein',ncs$aachange,"in.group","GENE.vars Low Impact Mutations","Low Impact");
p9 <- graph_variant_list(GENE.vars.mr,'PredictProtein','CADD.phred',ncs$aachange,"in.group","Low Impact Mutations","Low Impact");
####

## PDF output 2
pdf(file="outputs/GENE.vars/GENE.varscontrols3.pdf", height = 12, width = 17)
grid.arrange(tableGrob(lit), top = "GENE.vars Literature Mutations" );p01;p02;p03
grid.arrange(tableGrob(scs), top = "GENE.vars Synonymous Mutations" );p1;p2;p3
grid.arrange(tableGrob(pcs), top = "GENE.vars High Impact Mutations" );p4;p5;p6
grid.arrange(tableGrob(ncs), top = "GENE.vars Low Impact Mutations");p7;p8;p9;p9
#p10
dev.off()


#### PDF output
# pdf(file="outputs/GENE.vars/GENE.varscontrols.pdf", height = 12, width = 17)
# grid.arrange(tableGrob(scs), top = "GENE.vars Synonymous Controls" );p11;p12
# grid.arrange(tableGrob(pcs), top = "GENE.vars Positive Controls" );p1;p2;p3
# grid.arrange(tableGrob(ncs), top = "GENE.vars Negative Controls");p4;p5;p6;p9
# #p10
# dev.off()
####

#### What happens when we compare same aa changes (different nuc change)?
GENE.vars.mr.nc %>% filter(duplicated(aachange,fromLast=FALSE) | duplicated(aachange,fromLast=TRUE)) %>%
  arrange(aapos.x)-> sameaa
ggplot(sameaa, aes(x = aapos.x, y = CADD.phred)) + geom_point(aes(colour=PredictProtein)) + 
  xlab("AA Pos") + ylab("CADD Phred") + ggtitle("GENE.vars Same aa vars") + 
  scale_color_gradient2(low = "green", mid = "yellow", high = "red" )
####

# cor(GENE.vars.mr$LJB_MutationTaster, GENE.vars.mr$CADD.phred)

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




