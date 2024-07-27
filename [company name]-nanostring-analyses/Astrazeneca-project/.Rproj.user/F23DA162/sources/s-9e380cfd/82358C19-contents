#Identify th

#load in libraries
Packages <- c('tidyverse','ggplot2','openxlsx','RColorBrewer')
lapply(Packages, library, character.only = TRUE)

#library
library('tidyverse')
library('ggplot2')
library('openxlsx')
library('RColorBrewer')
library(rstatix)
library(ggpubr)
library(ggpmisc)
library(DBI)
library(dplyr)
library(dbplyr)
library(odbc)
library(readr)

#function
read_all_sheets = function(xlsxFile, ...) {
  sheet_names = openxlsx::getSheetNames(xlsxFile)
  sheet_list = as.list(rep(NA, length(sheet_names)))
  names(sheet_list) = sheet_names
  for (sn in sheet_names) {
    sheet_list[[sn]] = openxlsx::read.xlsx(xlsxFile, sheet=sn, ...)
  }
  return(sheet_list)
}

#cell type
cell_type_relative <- read.csv2("/home/bionerd/Dropbox/Company/Xsphera/@Data/Nanostring/Pooled/Full_Probe_Annotation 2021-03-19 09-56/results/cell types/cell type scores - relative.csv",
                       sep=",")
#TIL
cell_type_raw <- read.csv2("/home/bionerd/Dropbox/Company/Xsphera/@Data/Nanostring/Pooled/Full_Probe_Annotation 2021-03-19 09-56/results/cell types/cell type scores - raw.csv",
                       sep=",")

#pathway
pathway <- read.csv2("/home/bionerd/Dropbox/Company/Xsphera/@Data/Nanostring/Pooled/Full_Probe_Annotation 2021-03-19 09-56/results/pathway scoring/signature scores.csv",
                           sep=",")

#genes
expression <- read_all_sheets("/home/bionerd/Dropbox/Company/Xsphera/@Data/Nanostring/Pooled/Full_Probe_Annotation 2021-03-19 09-56/results/Normalization/mRNA_normalized_data.xlsx")

#metadata
metadata <- read.csv2('/home/bionerd/Dropbox/Company/Xsphera/@Data/Nanostring/Xsphera Server Nanostring Data/2021.03.11 Nanostring_Metadata.csv',
                      sep=",",fileEncoding="latin1")

#find missing sample names
#cell_type_raw %>% filter(!(Sample.Name %in% metadata$Sample.Name)) %>% select(Sample.Name)

#verify all samples can be found from charts
cell_type_raw %>% filter(Sample.Name %in% metadata$Sample.Name) %>% nrow()
pathway %>% filter(Sample.Name %in% metadata$Sample.Name) %>% nrow()
cell_type_relative %>% filter(Sample.Name %in% metadata$Sample.Name) %>% nrow()

#gather appropiate data & metadata for cell type
tmp <- cell_type_relative %>% select(Sample.Name, Total.TILs)
cell_type <- merge(cell_type_raw, tmp)
tmp <- metadata %>% select(Sample.Name,Disease,Drug,Day,Dosage,Patient.ID,pooled)
cell_type <- merge(cell_type, tmp)
cell_type <- cell_type %>% gather(Cell_Type,Value,2:15)

#gather appropiate data & metadaa for pathway
tmp <- metadata %>% select(Sample.Name,Disease,Drug,Day,Dosage,Patient.ID,pooled)
pathway <- merge(pathway, tmp)
pathway <- pathway %>% gather(pathway_name,Value,2:26)

#fix data types
cell_type$Value <- as.numeric(cell_type$Value)
pathway$Value <- as.numeric(pathway$Value)

#cell type distribution across all

#save plot

pdf(file="cell_type.pdf")
ggplot(cell_type, aes(x = Cell_Type, y = Value, fill = pooled)) +    
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  coord_flip()+
  ggtitle("Cell Type Over All Diseases") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#pathway across all
pdf(file="pathway.pdf")
ggplot(pathway, aes(x = pathway_name, y = Value, fill = pooled)) +    
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  coord_flip()+
  ggtitle("Pathway Over All Diseases") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Plot for each Disease the cell type and pathway figures
path <- "/c/Users/Nathan/Dropbox/Company/Xsphera/@Analysis/Pooled vs Not"

#for each disease
for (i in unique(cell_type$Disease)){
  print(i)
  #split to for pathway
  tmp <- pathway %>% filter(Disease == i)

  #save plot
  file_name = paste(path,"/",i,"_disease_pathway.pdf",sep="")
  pdf(file=file_name)

  #pathway across all
  print(ggplot(tmp, aes(x = pathway_name, y = Value, fill = pooled)) +    
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    coord_flip() +
    ggtitle(i) +
    theme(plot.title = element_text(hjust = 0.5)))

  dev.off()
  
  #split to disease
  tmp <- cell_type %>% filter(Disease == i)
  
  #save plot
  file_name = paste(path,"/",i,"_disease_cell_type.pdf",sep="")
  pdf(file=file_name)
  
  #cell type distribution across all
  print(ggplot(tmp, aes(x = Cell_Type, y = Value, fill = pooled)) +    
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    coord_flip() +
    ggtitle(i) +
    theme(plot.title = element_text(hjust = 0.5)))
  
  dev.off()
}

#for each drug
for (i in unique(cell_type$Drug)){
  print(i)
  #split to for pathway
  tmp <- pathway %>% filter(Drug == i)
  
  #save plot
  file_name = paste(path,"/",i,"_drug_pathway.pdf",sep="")
  pdf(file=file_name)
  
  #pathway across all
  print(ggplot(tmp, aes(x = pathway_name, y = Value, fill = pooled)) +    
          geom_boxplot() + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
          coord_flip() +
          ggtitle(i) +
          theme(plot.title = element_text(hjust = 0.5)))
  
  dev.off()
  
  #split to disease
  tmp <- cell_type %>% filter(Drug == i)
  
  #save plot
  file_name = paste(path,"/",i,"_drug_cell_type.pdf",sep="")
  pdf(file=file_name)
  
  #cell type distribution across all
  print(ggplot(tmp, aes(x = Cell_Type, y = Value, fill = pooled)) +    
          geom_boxplot() + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
          coord_flip() +
          ggtitle(i) +
          theme(plot.title = element_text(hjust = 0.5)))
  
  dev.off()
}

#for each drug and disease combination
for (i in unique(cell_type$Disease)){
  for (n in unique(cell_type$Drug)){
    print(c(i,n))
    #split to for pathway
    tmp <- pathway %>% filter(Disease == i & Drug == n)
    
    #save plot
    file_name = paste(path,"/",i,"_",n,"_drug_pathway.pdf",sep="")
    pdf(file=file_name)
    
    #pathway across all
    print(ggplot(tmp, aes(x = pathway_name, y = Value, fill = pooled)) +    
            geom_boxplot() + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
            coord_flip() +
            ggtitle(paste(i,"_",n,"Pathway",sep="")) +
            theme(plot.title = element_text(hjust = 0.5)))
    
    dev.off()
    
    #split to for cell type
    tmp <- cell_type %>% filter(Disease == i & Drug == n)
    
    #save plot
    file_name = paste(path,"/",i,"_",n,"_drug_cell_type.pdf",sep="")
    pdf(file=file_name)
    
    #cell type distribution across all
    print(ggplot(tmp, aes(x = Cell_Type, y = Value, fill = pooled)) +    
            geom_boxplot() + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
            coord_flip() +
            ggtitle(paste(i,"_",n,"Cell Type",sep="")) +
            theme(plot.title = element_text(hjust = 0.5)))
    
    dev.off()
    
    
    
    
  }
}






