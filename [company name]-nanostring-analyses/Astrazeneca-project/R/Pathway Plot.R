#Identify th

#load in libraries
Packages <- c('tidyverse','ggplot2','openxlsx','RColorBrewer',
              'rstatix','ggpubr')
lapply(Packages, library, character.only = TRUE)

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

#pathway
pathway <- read.csv2("/home/bionerd/Dropbox/Company/Xsphera/@Data/Nanostring/Astrazeneca/Full_minusD0_QC 2021-04-28 16-17/results/pathway scoring/signature scores.csv",
                     sep=",")
names(pathway)[1] <- "File.Name"

#metadata
metadata <- read.csv2('/home/bionerd/Dropbox/Company/Xsphera/@Data/Nanostring/Xsphera Server Nanostring Data/2021.04.28 Nanostring_Metadata.csv',
                      sep=",",fileEncoding="latin1")

#find missing sample names
#cell_type_raw %>% filter(!(Sample.Name %in% metadata$Sample.Name)) %>% select(Sample.Name)

#verify all samples can be found from charts
pathway %>% filter(File.Name %in% metadata$File.Name) %>% nrow()

#gather appropiate data & metadaa for pathway
tmp <- metadata %>% select(File.Name,Disease,Drug,Day,Dosage,Patient.ID)
pathway <- merge(pathway, tmp)
pathway <- pathway %>% gather(pathway_name,Value,2:26)

#fix data types
pathway$Value <- as.numeric(pathway$Value)

#pathway across all
pdf(file="pathway.pdf")
ggplot(pathway, aes(x = pathway_name, y = Value, fill = Drug)) +    
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  coord_flip()+
  ggtitle("Pathway Over All Drugs") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Plot for each Disease the cell type and pathway figures
path <- "/c/Users/Nathan/Dropbox/Company/Xsphera/@Analysis/Astrazeneca"

#for each Patient ID
for (i in unique(pathway$Patient.ID)){
  print(i)
  #split to for pathway
  tmp <- pathway %>% filter(Patient.ID == i)
  
  #save plot
  file_name = paste(path,"/",i,"_drug_pathway.pdf",sep="")
  pdf(file=file_name)
  
  #pathway across all
  print(ggplot(tmp, aes(x = pathway_name, y = Value, fill = Drug)) +    
          geom_boxplot() + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
          coord_flip() +
          ggtitle(i) +
          theme(plot.title = element_text(hjust = 0.5)))
  
  dev.off()
}
