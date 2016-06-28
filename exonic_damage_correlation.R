# SETUP -------------------------------------------------------------------
library(dplyr)
library(ssh.utils)
library(stringr)
library("biomaRt")
library(seqinr)
library(ggplot2)

source("/home/bcallaghan/NateDBCopy/functions.R")
source("/home/bcallaghan/NateDBCopy/exonic_damage_correlation_functions.R")

setwd("/home/bcallaghan/NateDBCopy/20kcorr/")
INPUTDIR <- "/home/bcallaghan/NateDBCopy/20kcorr/"


# BIOMART SETUP -----------------------------------------------------------

# GENE <- "GRIN2B"

cat("Starting analysis for", GENE,"\n")

cat("Loading Biomart...","\n")
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host="grch37.ensembl.org",path="/biomart/martservice")
BMattr <- listAttributes(mart)
# predictProteinMap <- read.table("/misc/pipeline42/ppdatabases/snap2results/UP000005640_9606/fastamap")
predictProteinMap <- read.table("fastamap")

BMtranscript.sort <- bioGetTranscript(GENE)
TRANSCRIPT <- BMtranscript.sort$refseq_mrna
cat("Done","\n")

# AMINO -------------------------------------------------------------------

cat("Loading snap2 results...","\n")
snapfile <- predictProteinMap[which(predictProteinMap$V2 == GENE),1]
fasta <- read.table(paste0("/misc/pipeline42/ppdatabases/snap2results/UP000005640_9606/",snapfile), skip=1)
snap2path <- paste0('/misc/pipeline42/ppdatabases/snap2results/UP000005640_9606/', snapfile, '.snap2.parsed')
snap2results  <- read.table(snap2path)
# compare_fastas(translated_cdna=cdna.translated, pp_file=GENE.PP, bm_fasta=FASTA)
cat("Done","\n")

# DNA ---------------------------------------------------------------------

cat("Creating Annovar Input File...","\n")
anno_in_path <- paste0(INPUTDIR,"tmp/",GENE,"_anno_in")
anno_outpath <- paste0(INPUTDIR,"tmp/",GENE)
anno_out_path <- paste0(INPUTDIR, "tmp/", GENE, ".hg19_multianno.csv")

vcf <- retrieve_exonic_vcf(GENE) # Retrieve vcf input for annovar for any gene
write.table(vcf, file = anno_in_path, quote = FALSE,
            sep="\t", row.names = FALSE, col.names = FALSE)

cat("Annovar input written to", anno_in_path,"\n")
cat("Done","\n")

cat("Running Annovar...","\n")
if(!file.exists(anno_out_path)){
  annocmd <- paste0("perl /space/bin/annovar/table_annovar.pl ", anno_in_path, " /space/bin/annovar/humandb/ -buildver hg19 -out ", anno_outpath, " -remove -protocol refGene,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp135,ljb_all,exac03,cadd -otherinfo -operation g,r,f,f,f,f,f,f -nastring . -csvout")
  annocmd
  cmd.out <- run.remote(cmd=annocmd , remote= "apu")
}
cat("Done","\n")

# ANNOTATED ---------------------------------------------------------------

multianno <- read.csv(anno_out_path)
head(multianno)


# MERGE -------------------------------------------------------------------

cat("Merging damage files...","\n")
multianno.2 <- merge_multianno_and_snap2(multianno, snap2results)
cat("Done","\n")

# CORRELATIONS ------------------------------------------------------------

cat("Calculating Damage Correlations...","\n")
genecor <- paste0(GENE,"\t",cor(multianno.2$CADD.phred, multianno.2$snap2, method="spearman"))
write(genecor,"outputs/genecorrelations", append=TRUE)
cat("Done","\n")

# OUTPUTS -----------------------------------------------------------------
cat("Outputting plots...","\n")
save_plots(multianno.2,"correlation")
save_plots(multianno.2,"snap2")
save_plots(multianno.2,"cadd")
# pngoutputpath <- paste0(INPUTDIR, "outputs/", GENE, "_CScor.png")
# p1 <- ggplot(multianno.2,aes(x= CADD.phred, y = snap2)) + geom_point()
# ggsave(pngoutputpath, p1)
# png(pngoutputpath,width = 480, height = 480)
# p1
# dev.off()
cat("Done","\n")

cat("Finished Analysis for:", GENE,"\n")
cat(paste0("Outputs in ", INPUTDIR, "tmp/","\n"))

setwd("/home/bcallaghan/NateDBCopy/")