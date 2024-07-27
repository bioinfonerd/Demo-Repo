
#TIS Analysis
#Tumor Inflammatory Signature on Astrazeneca


library(stringr)


#data
TIS.signature <- c("CCL5-mRNA", "CD8A-mRNA","STAT1-mRNA","CXCL9-mRNA","CD27-mRNA","CXCR6-mRNA","IDO1-mRNA",
                   "TIGIT-mRNA","LAG3-mRNA","CD276-mRNA","PDCD1LG2-mRNA","CD274-mRNA","HLA-E-mRNA","HLA-DQA1-mRNA",
                   "PSMB10-mRNA","HLA-DRB1-mRNA","CMKLR1-mRNA","NKG7-mRNA")

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
metadata <- read.csv2('/home/bionerd/Dropbox/Company/Xsphera/@Data/Nanostring/Xsphera Server Nanostring Data/2021.04.28 Nanostring_Metadata.csv',
                      sep=",",fileEncoding="latin1")

#calculate TIS Score for each sample with or without Durva
results <- expression$Log2 %>% select(Sample.Name,TIS.signature)
results$TIS_Sig_Sum <- rowSums(results[ , c(2,length(TIS.signature))], na.rm=TRUE)
results$TIS_Sig_Mean <- rowMeans(results[ , c(2,length(TIS.signature))], na.rm=TRUE)

#add meta information
tmp <- metadata %>% select(Patient.ID,Sample.Name,Drug,Dosage)
output <- merge(results,tmp,by.x="Sample.Name",by.y="Sample.Name")

#Plot 1 TIS Score - Overall (SUM)
pdf(file="./images/TIS.Score.VSDrug.Sum.boxplot.bin.pdf")
ggplot(output, aes(x=Drug, y=TIS_Sig_Sum)) +
  geom_boxplot()+
  labs(title = "TIS.Score By Drug",
       x = "Drug",
       y = "Sum of TIS Signature Genes")
dev.off()

# Plot 1b TIS Score - Overall (MEAN)
pdf(file="./images/TIS.Score.VSDrug.Mean.boxplot.bin.pdf")
ggplot(output, aes(x=Drug, y=TIS_Sig_Mean)) +
  geom_boxplot()+
  labs(title = "TIS.Score By Drug",
       x = "Drug",
       y = "Mean of TIS Signature Genes")
dev.off()

#Plot 2 - concentration
tmp <- output %>% filter(Drug=="IgG" | Drug == "Durva")
tmp <- tmp %>% filter(Patient.ID=="P0113")

pdf(file="./images/TIS.Score.VSDrug.concentration.Mean.boxplot.bin.pdf")
ggplot(tmp, aes(x=Dosage, y=TIS_Sig_Mean)) +
  geom_boxplot(aes(colour = Drug))+
  labs(title = "TIS.Score By Drug",
       x = "Concentration",
       y = "Mean of TIS Signature Genes")
dev.off()

#Plot 3 Nomacan 30 uG vs IgG
expression <- read_all_sheets('/c/Users/Nathan/Dropbox/Company/Xsphera/@Data/Nanostring/Nomacan_Meta_Analysis/Full_IgG_Probe 2021-03-05 12-11/results/Normalization/Full-mRNA_normalized_data.xlsx')
results <- expression$Sheet1 %>% select(File.Name,TIS.signature)
results[ , c(2,length(TIS.signature))] <- apply(results[ , c(2,length(TIS.signature))], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))
results$TIS_Sig_Sum <- rowSums(results[ , c(2,length(TIS.signature))], na.rm=TRUE)
results$TIS_Sig_Mean <- rowMeans(results[ , c(2,length(TIS.signature))], na.rm=TRUE)
tmp <- metadata %>% select(Patient.ID,File.Name,Drug,Dosage)
output <- merge(results,tmp,by.x="File.Name",by.y="File.Name")


tmp <- output %>% filter(Drug=="IgG" | Drug == "NMC2" | Drug == "NMC1")
tmp <- tmp %>% filter(Patient.ID=="P0113")

pdf(file="./images/TIS.Score.VSNomacanDrug.Mean.boxplot.bin.pdf")
ggplot(tmp, aes(x=Drug, y=TIS_Sig_Mean)) +
  geom_boxplot(aes(colour = Drug))+
  labs(title = "TIS.Score By Drug",
       x = "Concentration",
       y = "Mean of TIS Signature Genes")
dev.off()


pdf(file="./images/TIS.Score.VSNomacanDrug.Concentration.Mean.boxplot.bin.pdf")
ggplot(tmp, aes(x=Dosage, y=TIS_Sig_Mean)) +
  geom_boxplot(aes(colour = Drug))+
  labs(title = "TIS.Score By Drug",
       x = "Concentration",
       y = "Mean of TIS Signature Genes")
dev.off()



