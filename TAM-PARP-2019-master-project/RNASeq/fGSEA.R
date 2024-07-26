## ----echo=F, message=F---------------------------------------------------
library(fgsea)
library(data.table)
library(ggplot2)
library(sleuth)
library(GSA)
library(dplyr)

# Location ----------------------------------------------------------------
setwd('current_directory')
source("gene_ids_entrez_conversion.R")

#gmt <- list.files("./GSEA", pattern = "\\.entrez.gmt$")
gmt <- list.files("./GSEA", pattern = "\\.symbols.gmt$")

# Analysis ---------------------------------------------------------------
files <- list.files(pattern = "\\.sleuthobj$")

i=5

sleuth<-sleuth_load(files[i])

#results table
table <- sleuth_results(sleuth,test = 'conditionOLA',
                        which_model = 'full',test_type = 'wt', show_all = TRUE, pval_aggregate = FALSE)

#table <- sleuth_results(sleuth,test = 'DayD9',
#                        which_model = 'full',test_type = 'wt', show_all = FALSE, pval_aggregate = FALSE)

table <- sleuth_results(sleuth, 'reduced:full', test_type = 'lrt',show_all = FALSE, pval_aggregate = FALSE)

matrix <- sleuth_to_matrix(sleuth,"obs_norm","est_counts")

s2c <- read.table(file.path("../", "Sample_Sheet.csv"), header = TRUE,sep=',', stringsAsFactors=FALSE) #metadata


CASP3
CASP7
CASP9
MRE11A
PARG
POLB
PARP1
PCNA
RFC4
RPA2
SNAI1




#grab 
rowMeans(matrix,)

#filter out duplicate genes (transcripts) & keep lowest p-value or [TODO] highest gene expression
table <- table %>%
  group_by(ext_gene) %>%
  slice(which.min(pval))

#convert gene ids to entrez
#id_conversion<-convert_symbol_entrez_id(table$ext_gene) #run conversion
#table <- table %>% #filter by symbols able to be converted
#  filter(ext_gene %in% id_conversion$query)
#table['entrez_id']<-id_conversion$entrez_id

#named vector of ranked gene by p-value (entrez)
#rank_gene_pval <- as.vector(t(table$b)) #fill vector with numeric
#names(rank_gene_pval) <- c(table$entrez_id) #id each element


#named vector of ranked gene by p-value (entrez)
rank_gene_pval <- as.vector(t(table$b)) #fill vector with numeric
names(rank_gene_pval) <- c(table$ext_gene) #id each element

#for D5,D6,D9 analysis since no b value
#rank_gene_pval <- as.vector(t(table$pval)) #fill vector with numeric
#names(rank_gene_pval) <- c(table$entrez_id) #id each element


#raw expression values


for (n in 1:length(gmt)){

gsea<-gmtPathways(paste(c('./GSEA/',gmt[n]),collapse = ''))
  
#gene_list<-GSA.read.gmt(paste(c('./GSEA/',gmt[n]),collapse = ''))
#gsea <- as.vector(gene_list$genesets)
#names(gsea) <- c(gene_list$geneset.names) 

#run analysis
output <- fgsea(pathways = gsea,
                       stats = rank_gene_pval,
                       minSize=15,
                       maxSize=500,
                       nperm=10000)

## ---- fig.width=10, fig.height=10, fig.retina=2----------------------------
topPathwaysUp <- output[ES > 0][head(order(pval), n=8), pathway]
topPathwaysDown <- output[ES < 0][head(order(pval), n=8), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

jpeg(paste(c(files[i],'.',gmt[n],'.topenriched.jpg'),collapse = ''))
plotGseaTable(gsea[topPathways], rank_gene_pval, output,
              colwidths=c(5, 3,0.8, 1.2, 1.2),gseaParam = 0.5)
dev.off()

#write data
fwrite(output, file=paste(c(files[i],'.',gmt[n],'.tsv'),collapse = ''), sep="\t", sep2=c("", " ", ""))
}
