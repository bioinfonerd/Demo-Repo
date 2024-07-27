#Goal To Repeat the complete Nanostring Nsolver Pipeline
#Addition of different methods to pipeline
#Nathan T. Johnson
#Edited: 2021.01.13 

# Libraries (Install / Load) ----------------------------------------------
#library("devtools")
#install_github("dbdimitrov/BingleSeq")
#library(BingleSeq)
#startBingleSeq()


#Verify Package Manager Installed, if not install
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", dependencies = TRUE)


#load in libraries
packages <- c('NanoStringDiff','NanoStringNorm','tidyverse','ggplot2','NACHO',
              'openxlsx','RColorBrewer','Hmisc','hrbrthemes','viridis',
              'knitr','PCAtools','biomaRt','corrplot')
#remotes::install_github("mcanouil/NACHO")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      BiocManager::install(x, ask = FALSE, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)




# Functions ---------------------------------------------------------------

#read in all sheets from excel
read_all_sheets = function(xlsxFile, ...) {
  sheet_names = openxlsx::getSheetNames(xlsxFile)
  sheet_list = as.list(rep(NA, length(sheet_names)))
  names(sheet_list) = sheet_names
  for (sn in sheet_names) {
    sheet_list[[sn]] = openxlsx::read.xlsx(xlsxFile, sheet=sn, ...)
  }
  return(sheet_list)
}

#function for saving
save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#loading in column information with integars adds an X (fix)
Correct_Colnames <- function(df) {
  
  #if X added to column, should not be empty
  delete.columns <- grep("(^X)", colnames(df), perl=T)
  
  #if not empty, replace 
  if (length(delete.columns) > 0) {
    colnames(df) <- gsub("(^X)", "",  colnames(df))    
  }
  
  return(df)
}



# RCC Data (Skipped For Now) ---------------------------------------------------

# #load in meta data 
# meta <- read.table('/home/bionerd/Dropbox/Company/Xsphera/@Data/Nanostring/Nanostring_Metadata.csv',
#                    sep=",",header=TRUE, stringsAsFactors=FALSE, fileEncoding="latin1")
# 
# #select meta data for only samples to be loaded
# files <- list.files('/home/bionerd/Dropbox/Company/Xsphera/@Data/Nanostring/Xsphera Server Nanostring Data/2021.01.31-Complete_RCC')
# 
# #filter
# meta <- meta %>% filter(File.Name %in% files)
# 
# #load in RCC data
# rcc_data <- load_rcc('/home/bionerd/Dropbox/Company/Xsphera/@Data/Nanostring/Xsphera Server Nanostring Data/2021.01.31-Complete_RCC',
#          ssheet_csv = meta,id_colname="File.Name")






# Load Data: Starting from nSolver Normalized Data -----------------------------

#base path
base='/home/bionerd/Dropbox/Company/Xsphera'

#set working directory
setwd(file.path(base,'@Analysis/Nomacan-Metaanalysis'))

#Pull in Nanostring metadata and expression data
Metadata <- read.csv(file.path(base,'@Data/Nanostring/Xsphera Server Nanostring Data/2021.01.14 Nanostring_Metadata=Nomacan Focused.csv'))
Expression_Data <- read.csv(file.path(base,'@Data/Nanostring/Nomacan_Meta_Analysis/Nomacan Meta Analysis_NormalizedData.csv'))
Genes <- read.csv(file.path(base,'@Data/Nanostring/Nanostring_LBL-10540-02_nCounter_NHP_Immunology_V2_Panel_Gene_List.csv'))
Annotations <- read.csv(file.path(base,'@Data/Nanostring/Nanostring_LBL-10540-02_nCounter_NHP_Immunology_V2_Panel_Gene_List_Annotations.csv'))
nanostring_cell_type_scores <- read.csv('/home/bionerd/Dropbox/Company/Xsphera/@Data/Nanostring/Media/New Analysis 2021-02-26 11-46/results/cell types/cell type scores - raw.csv')

#R puts an X before any column name that starts with integar in dataframe, fix
Expression_Data <- Correct_Colnames(Expression_Data)

#Flow Data (need % cells tab, sheet 1 = raw data)
flow_data <- read.xlsx(file.path(base,'@Data/Flow Data/20210218 baseline P0003-113-Nathan Modified.xlsx'),
                       sheet = 2)

# #pull in Pathway data
# full.pathway <- read.csv(file.path(base,'/@Data/Nanostring/Nomacan_Meta_Analysis/Full_IgG_Probe 2021-03-05 12-11/results/pathway scoring/signature scores.csv'))
# RCC.pathway <- read.csv(file.path(base,'/@Data/Nanostring/Nomacan_Meta_Analysis/RCC_IgG_probe 2021-03-05 12-06/results/pathway scoring/signature scores.csv'))
# lung.pathway <- read.csv(file.path(base,'/@Data/Nanostring/Nomacan_Meta_Analysis/lung_IgG_probe 2021-03-05 12-05/results/pathway scoring/signature scores.csv'))
# Melanoma.pathway <- read.csv(file.path(base,'/@Data/Nanostring/Nomacan_Meta_Analysis/Melanoma_IgG_Probe 2021-03-05 11-59/results/pathway scoring/signature scores.csv'))
# 
# #pull in GSEA data
# full.GSEA.directed <- read.csv(file.path(base,'/@Data/Nanostring/Nomacan_Meta_Analysis/Full_IgG_Probe 2021-03-05 12-11/results/pathway scoring/signature scores.csv'))
# RCC.GSEA.directed <- read.csv(file.path(base,'/@Data/Nanostring/Nomacan_Meta_Analysis/RCC_IgG_probe 2021-03-05 12-06/results/pathway scoring/signature scores.csv'))
# lung.GSEA.directed <- read.csv(file.path(base,'/@Data/Nanostring/Nomacan_Meta_Analysis/lung_IgG_probe 2021-03-05 12-05/results/pathway scoring/signature scores.csv'))
# Melanoma.GSEA.directed <- read.csv(file.path(base,'/@Data/Nanostring/Nomacan_Meta_Analysis/Melanoma_IgG_Probe 2021-03-05 11-59/results/pathway scoring/signature scores.csv'))
# 




# Flow Data (Correlate Cell Type Flow Data to D1 Media)------------------------------

#gather data into column
flow_data <- flow_data %>% gather("cell_type", "percent", -Sample, -fraction)
#replace NA with 0
flow_data[is.na(flow_data)] = 0

#plot violin plot of data
p <- ggplot(flow_data, aes(x=cell_type, y=percent)) + 
  geom_violin(trim=FALSE)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_boxplot(width=0.1, fill="white")
p

ggsave('Flow Data.svg',
  plot = p)

#calculate the relative % of cell type ()
#how to calculate total RNA amount in Nanostring
#total sum of counts of genes / total sum of counts of genes related to cell types
# 


#format data, select only media patient data and gene name
Expression_Data_flow <- Expression_Data %>% dplyr::select(-c("Annotation",
                                        "Accession..",
                                        "NS.Probe.ID",
                                        "Class.Name",
                                        "Analyte.Type",
                                        "..Samples.above.Threshold"),)

#remove positive and negative control
Expression_Data_flow <- Expression_Data_flow %>% filter(!grepl("NEG",Probe.Name) & !grepl("POS",Probe.Name))

#rotate expression data so genes are columns and clean it up
Expression_Data_flow <- as_tibble(cbind(nms = names(Expression_Data_flow), t(Expression_Data_flow)))
colnames(Expression_Data_flow) <-Expression_Data_flow[1,]
Expression_Data_flow <- Expression_Data_flow[-1, ]
#rename column
colnames(Expression_Data_flow)[colnames(Expression_Data_flow) == 'Probe.Name'] <- 'Sample.Name'

#filter to media only data
flow_sample_comparison <- Metadata %>% filter(Drug == 'Media' & Day == 'D0')
Expression_Data_flow <- Expression_Data_flow %>% filter(Sample.Name %in% flow_sample_comparison$File.Name)

#turn into numeric and calculate sum across each row that is numeric (in order to exclude sample name column)
Expression_Data_flow[,-1] <- lapply(Expression_Data_flow[,-1], function(x) as.numeric(x))
Expression_Data_flow <- Expression_Data_flow %>% 
                            rowwise(Sample.Name) %>% 
                            mutate(row.sum = sum(c_across(where(is.numeric))))

#Cell Types Available
Cell.Types <- unique(Annotations$Cell.Type)

#make tibble for sum cell type
# celltype_rna_expression_sum <- tibble("Sample.Name" = Expression_Data_flow$Sample.Name)
# for (i in 2:length(Cell.Types)){
#   selected_genes <- Annotations %>% filter(Cell.Type == Cell.Types[i])  
#   tmp <- Expression_Data_flow %>%
#     rowwise(Sample.Name) %>%
#     dplyr::select(matches(selected_genes$Gene)) %>%
#     mutate(!!Cell.Types[i] := sum(c_across(where(is.numeric))))
#   celltype_rna_expression_sum[paste(Cell.Types[i],"nanostring")] <- tmp[Cell.Types[i]]
# }
# 
# #make tibble for average cell type
# celltype_rna_expression_avg <- tibble("Sample.Name" = Expression_Data_flow$Sample.Name)
# for (i in 2:length(Cell.Types)){
#   selected_genes <- Annotations %>% filter(Cell.Type == Cell.Types[i])  
#   tmp <- Expression_Data_flow %>%
#     rowwise(Sample.Name) %>%
#     dplyr::select(matches(selected_genes$Gene)) %>%
#     mutate(!!Cell.Types[i] := mean(c_across(where(is.numeric))))
#   celltype_rna_expression_avg[paste(Cell.Types[i],"nanostring")] <- tmp[Cell.Types[i]]
# }
# 
# #make tibble for max cell type
# celltype_rna_expression_max <- tibble("Sample.Name" = Expression_Data_flow$Sample.Name)
# for (i in 2:length(Cell.Types)){
#   selected_genes <- Annotations %>% filter(Cell.Type == Cell.Types[i])  
#   tmp <- Expression_Data_flow %>%
#     rowwise(Sample.Name) %>%
#     dplyr::select(matches(selected_genes$Gene)) %>%
#     mutate(!!Cell.Types[i] := max(c_across(where(is.numeric))))
#   celltype_rna_expression_max[paste(Cell.Types[i],"nanostring")] <- tmp[Cell.Types[i]]
# }


#compile results from flow with nanostring (1 flow represention per sample, but multiple for nanostring)
names(nanostring_cell_type_scores)[2:length(nanostring_cell_type_scores)]

#
nanostring_cell_type_scores$CD45
nanostring_cell_type_scores$Sample.Name

#filter for cell type
#filter for sample name

sample_names<-c("P0092","P0092","P0092","P0071","P0071","P0071","P0071",
"P0071","P0071","P0083","P0083", "P0083", "P0094", "P0094",
"P0094", "P0090","P0090", "P0090", "P0091", "P0091",       
"P0091", "P0096", "P0096", "P0096", "P0096", "P0096","P0096", "P0103", "P0103", "P0103")   

cd45<-1:nrow(nanostring_cell_type_scores)
for (i in 1:nrow(nanostring_cell_type_scores)){
  cd45[i] <-   flow_data %>%
    filter(Sample == sample_names[i] & fraction == "S1" & cell_type == "CD45+") %>%
    dplyr::select(percent)
}

tmp <- flow_data %>%
  filter(Sample == sample_names[i] & fraction == "S1" & cell_type == "CD45+") %>%
  dplyr::select(percent)

nanostring_cell_type_scores['cd45_flow'] <- cd45
  
cbind(nanostring_cell_type_scores, cd45_flow=cd45)


cbind(nanostring_cell_type_scores, cd45)

vec <- c(3, 2, 3, 2, 3) 
#plot cell type data
library("ggpubr")
ggscatter(x = "mpg", y = cd45, 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Miles/(US) gallon", ylab = "Weight (1000 lbs)")





#for each patient, capture flow data and media
#store data as correlation value encoded for each cell type and patient

#capture genes for cell type
#calculate % from RNA expression?







# PCA ---------------------------------------------------------------------

# Can we group the samples?

#expression data formatting
#remove gene controls from expression data
#PCA_Data <- Expression_Data %>% filter(!grepl("NEG",Gene_Name) & !grepl("POS",Gene_Name))
#row.names(PCA_Data) <- PCA_Data$Gene_Name

#drop unrelevant columns
Genes <- Expression_Data$Probe.Name
PCA_Data <- Expression_Data %>% dplyr::select(-c("Probe.Name",
                                                 "Annotation",
                                                 "Accession..",
                                                 "NS.Probe.ID",
                                                 "Class.Name",
                                                 "Analyte.Type",
                                                 "..Samples.above.Threshold"),)

PCA_Data <- as.matrix((PCA_Data))

#metadata formatting
PCA_Metadata <- Metadata %>% dplyr::select(c("File.Name","Patient.ID",
                                             "Disease","Drug","Day"))
row.names(PCA_Metadata) <- Metadata$File.Name

#Conduct PCA
p <- pca(PCA_Data, metadata = PCA_Metadata, removeVar = 0.1)


identical(colnames(PCA_Data), rownames(PCA_Metadata))
if (!identical(rownames(PCA_Data), colnames(PCA_Metadata))) {
  stop("'colnames(mat)' is not identical to 'rownames(metadata)'")
}


if (!identical(colnames(PCA_Data), rownames(PCA_Metadata))) {
  stop("'colnames(mat)' is not identical to 'rownames(metadata)'")
}

colnames(PCA_Data)[1]
rownames(PCA_Metadata)[1]
setdiff(rownames(PCA_Metadata),colnames(PCA_Data))

(new <- rownames(PCA_Metadata)[!rownames(PCA_Metadata) %in% colnames(PCA_Data)])
(new <- colnames(PCA_Metadata)[!colnames(PCA_Metadata) %in% rownames(PCA_Data)])

#Plot screen plot
output <- screeplot(p,returnPlot = TRUE)
ggplot2::ggsave('C:/Users/Nathan/Dropbox/Company/Xsphera/@Data/Nanostring/Analysis/PCAPlot.png')

#Plot PC1 vs PC2
biplot(p,    
       colby = 'Tumor',
       hline = 0, vline = 0,
       legendPosition = 'right', legendLabSize = 16, legendIconSize = 8.0,
       shape = 'Treatment',
       drawConnectors = FALSE,
       title = 'PC1 vs PC2',
       subtitle = 'PC1 versus PC2')
ggplot2::ggsave('C:/Users/Nathan/Dropbox/Company/Xsphera/@Data/Nanostring/Analysis/Pc1vsPC2.png')


#Plot each component
pairsplot(p,
          triangle = FALSE,
          hline = 0, vline = 0,
          pointSize = 0.8,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'Tumor',
          title = 'Pairs plot', plotaxes = TRUE,
          margingaps = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))
ggplot2::ggsave('C:/Users/Nathan/Dropbox/Company/Xsphera/@Data/Nanostring/Analysis/Principle Component Comparison.png')

#All Plots
eigencorplot(p,
             metavars = colnames(PCA_Metadata),
             col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
             cexCorval = 0.7,
             colCorval = 'white',
             fontCorval = 2,
             posLab = 'bottomleft',
             rotLabX = 45,
             posColKey = 'top',
             cexLabColKey = 1.5,
             scale = TRUE,
             main = 'Principle Components Correlation',
             colFrame = 'white',
             plotRsquared = FALSE)
ggplot2::ggsave('C:/Users/Nathan/Dropbox/Company/Xsphera/@Data/Nanostring/Analysis/Principle Component MetaData Correlation.png')

#Which Principle Component is significantly Correlated 
eigencorplot(p,
             metavars = colnames(PCA_Metadata),
             col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
             cexCorval = 0.7,
             colCorval = 'white',
             fontCorval = 2,
             posLab = 'bottomleft',
             rotLabX = 45,
             posColKey = 'top',
             cexLabColKey = 1.5,
             scale = TRUE,
             main = bquote(Principal ~ component ~ Pearson ~ r^2 ~ clinical ~ correlates),
             colFrame = 'white',
             plotRsquared = TRUE,
             corFUN = 'pearson',
             corUSE = 'pairwise.complete.obs',
             signifSymbols = c('****', '***', '**', '*', ''),
             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1))
ggplot2::ggsave('C:/Users/Nathan/Dropbox/Company/Xsphera/@Data/Nanostring/Analysis/Signficiant Principle Component MetaData Correlation.png')


#What genes are driving the principle component
plotloadings(p,
             rangeRetain = 0.01,
             labSize = 3.0,
             title = 'What genes are driving the principle component',
             subtitle = 'PC1, PC2, PC3, PC4, PC5',
             caption = 'Top 1% variables',
             shape = 24,
             col = c('limegreen', 'black', 'red3'),
             drawConnectors = TRUE)
ggplot2::ggsave('C:/Users/Nathan/Dropbox/Company/Xsphera/@Data/Nanostring/Analysis/Genes Driving Principle Compoments (1-5).png')

#What genes are driving the principle component
plotloadings(p,
             rangeRetain = 0.01,
             components = getComponents(p, c(1,3,4,9)),
             labSize = 3.0,
             title = 'What genes are driving the principle component',
             subtitle = 'Focus on Significant Components',
             shape = 24,
             col = c('limegreen', 'black', 'red3'),
             drawConnectors = TRUE)
ggplot2::ggsave('C:/Users/Nathan/Dropbox/Company/Xsphera/@Data/Nanostring/Analysis/Genes Driving Significant Principle Compoments.png')


#PC1 and PC3 seem to be the main areas of interest
#Plot PC1 vs PC2
biplot(p,    
       x = 'PC1',
       y = 'PC3',
       colby = 'Tumor',
       hline = 0, vline = 0,
       legendPosition = 'right', legendLabSize = 16, legendIconSize = 8.0,
       shape = 'Treatment',
       drawConnectors = FALSE,
       title = 'Component Comparison',
       subtitle = 'PC1 versus PC3')
ggplot2::ggsave('C:/Users/Nathan/Dropbox/Company/Xsphera/@Data/Nanostring/Analysis/Pc1vsPC3.png')


#ID gene functions [TODO]
# mart <- useMart('ENSEMBL_MART_ENSEMBL')
# mart <- useDataset('hsapiens_gene_ensembl', mart)
# 
#   getBM(mart = mart,
#     attributes = c('affy_hg_u133a', 'ensembl_gene_id',
#       'gene_biotype', 'external_gene_name'),
#     filter = 'affy_hg_u133a',
#     values = c('215281_x_at', '214464_at', '211122_s_at', '205225_at',
#       '202037_s_at', '204540_at', '215176_x_at', '205044_at', '208650_s_at',
#       '205380_at'),
#     uniqueRows = TRUE)





# Can we group the samples without treatment samples?

#metadata formatting & remove samples with treatment
sample_to_keep <- (Metadata %>% filter(Treatment == 'None'))$Sample.Name
PCA_Metadata <- Metadata %>% dplyr::select(-c("File.Name","Sample.Name","Cartridge.ID","Lane.Number","Import.Date")) %>% filter(Treatment == 'None')
row.names(PCA_Metadata) <- sample_to_keep

#expression data formatting
#remove gene controls from expression data 
PCA_Data <- Expression_Data %>% filter(!grepl("NEG",Gene_Name) & !grepl("POS",Gene_Name))
row.names(PCA_Data) <- PCA_Data$Gene_Name
#drop unrelevant columns & remove samples with treatment
PCA_Data <- PCA_Data %>% dplyr::select(-c("Accession","Sample.ID","Gene_Name"))
PCA_Data <- PCA_Data %>% dplyr::select(c(sample_to_keep))
PCA_Data <- as.matrix((PCA_Data))

#Conduct PCA
p <- pca(PCA_Data, metadata = PCA_Metadata, removeVar = 0.1)

#Plot scree plot
output <- screeplot(p,returnPlot = TRUE)
ggplot2::ggsave('C:/Users/Nathan/Dropbox/Company/Xsphera/@Data/Nanostring/Analysis/PCAPlot.png')

#Plot PC1 vs PC2
biplot(p,    
       colby = 'Tumor',
       hline = 0, vline = 0,
       legendPosition = 'right', legendLabSize = 16, legendIconSize = 8.0,
       shape = 'Treatment',
       drawConnectors = FALSE,
       title = 'PC1 vs PC2',
       subtitle = 'PC1 versus PC2')
ggplot2::ggsave('C:/Users/Nathan/Dropbox/Company/Xsphera/@Data/Nanostring/Analysis/Pc1vsPC2_notreatmentsamples.png')


#Plot each component
pairsplot(p,
          triangle = FALSE,
          hline = 0, vline = 0,
          pointSize = 0.8,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'Tumor',
          title = 'Pairs plot', plotaxes = TRUE,
          margingaps = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))
ggplot2::ggsave('C:/Users/Nathan/Dropbox/Company/Xsphera/@Data/Nanostring/Analysis/Principle Component Comparison_notreatmentsamples.png')

#All Plots
# eigencorplot(p,
#   metavars = colnames(PCA_Metadata),
#   col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
#   components = getComponents(p),
#   cexCorval = 0.7,
#   colCorval = 'white',
#   fontCorval = 2,
#   posLab = 'bottomleft',
#   rotLabX = 45,
#   posColKey = 'top',
#   cexLabColKey = 1.5,
#   scale = TRUE,
#   main = 'Principle Components Correlation',
#   colFrame = 'white',
#   plotRsquared = FALSE)
# ggplot2::ggsave('C:/Users/Nathan/Dropbox/Company/Xsphera/@Data/Nanostring/Analysis/Principle Component MetaData Correlation_notreatmentsamples.png')

# #Which Principle Component is significantly Correlated 
# eigencorplot(p,
#   metavars = colnames(PCA_Metadata),
#   col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
#   components = getComponents(p),
#   cexCorval = 0.7,
#   colCorval = 'white',
#   fontCorval = 2,
#   posLab = 'bottomleft',
#   rotLabX = 45,
#   posColKey = 'top',
#   cexLabColKey = 1.5,
#   scale = TRUE,
#   main = bquote(Principal ~ component ~ Pearson ~ r^2 ~ clinical ~ correlates),
#   colFrame = 'white',
#   plotRsquared = TRUE,
#   corFUN = 'pearson',
#   corUSE = 'pairwise.complete.obs',
#   signifSymbols = c('****', '***', '**', '*', ''),
#   signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1))
# ggplot2::ggsave('C:/Users/Nathan/Dropbox/Company/Xsphera/@Data/Nanostring/Analysis/Signficiant Principle Component MetaData Correlation_notreatmentsamples.png')


#What genes are driving the principle component
plotloadings(p,
             rangeRetain = 0.01,
             labSize = 3.0,
             title = 'What genes are driving the principle component',
             subtitle = 'PC1, PC2, PC3, PC4, PC5',
             caption = 'Top 1% variables',
             shape = 24,
             col = c('limegreen', 'black', 'red3'),
             drawConnectors = TRUE)
ggplot2::ggsave('C:/Users/Nathan/Dropbox/Company/Xsphera/@Data/Nanostring/Analysis/Genes Driving Principle Compoments (1-5)_notreatmentsamples.png')

#What genes are driving the principle component
# plotloadings(p,
#   rangeRetain = 0.01,
#   components = getComponents(p, c(1,3,4,9)),
#   labSize = 3.0,
#   title = 'What genes are driving the principle component',
#   subtitle = 'Focus on Significant Components',
#   shape = 24,
#   col = c('limegreen', 'black', 'red3'),
#   drawConnectors = TRUE)
# ggplot2::ggsave('C:/Users/Nathan/Dropbox/Company/Xsphera/@Data/Nanostring/Analysis/Genes Driving Significant Principle Compoments_notreatmentsamples.png')


#PC1 and PC3 seem to be the main areas of interest
#Plot PC1 vs PC2
#   biplot(p,    
#          x = 'PC1',
#          y = 'PC3',
#          colby = 'Tumor',
#          hline = 0, vline = 0,
#          legendPosition = 'right', legendLabSize = 16, legendIconSize = 8.0,
#          shape = 'Treatment',
#          drawConnectors = FALSE,
#          title = 'Component Comparison',
#          subtitle = 'PC1 versus PC3')
# ggplot2::ggsave('C:/Users/Nathan/Dropbox/Company/Xsphera/@Data/Nanostring/Analysis/Pc1vsPC3_notreatmentsamples.png')








# HeatMap & Clustermap ----------------------------------------------------



pheatmap::pheatmap(output,main=CD45_flow,annotation_row = samples)



# Pathway Scoring ---------------------------------------------------------


# Cell Type Scoring -------------------------------------------------------

#
#drop unrelevant columns
Genes <- Expression_Data$Probe.Name
data <- Expression_Data %>% dplyr::select(-c("Annotation",
                                                 "Accession..",
                                                 "NS.Probe.ID",
                                                 "Class.Name",
                                                 "Analyte.Type",
                                                 "..Samples.above.Threshold"),)

#rotate expression data so genes are columns and clean it up
data <- as_tibble(cbind(nms = names(data), t(data)))
colnames(data) <-data[1,]
data <- data[-1, ]

#meta information orientated to data
samples <- data.frame(data$Probe.Name)
colnames(samples) <- "File.Name"
samples <- samples %>% left_join(Metadata)
rownames(samples) <- samples$Sample.Name
samples <- samples %>% dplyr::select(c("File.Name","Patient.ID",
                                       "Disease","Drug","Day")) #remove unncessary colnames

#function for saving
save_pheatmap_png <- function(x, filename, width=3200, height=2000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#Cell Types Available
Cell.Types <- unique(Annotations$Cell.Type)

setwd('/c/Users/Nathan/Dropbox/Company/Xsphera/@Analysis/Nanostring')

for(i in 2:length(Cell.Types)){
  #Select genes available for Cell Type
  selected_genes <- Annotations %>% filter(Cell.Type == Cell.Types[i])
  
  if(nrow(selected_genes)>=2){
    #filter genes by cell type
    #not all genes are within dataset, exclude those
    output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
    if(length(colnames(output))>1){
      output <- sapply(output, as.numeric)
      rownames(output) <- rownames(samples)
      
      #plot
      my_heatmap <- pheatmap::pheatmap(output,
                                       main=Cell.Types[i],
                                       annotation_row = samples,
                                       scale="column")
      #save heatmap
      save_pheatmap_png(my_heatmap, paste(as.character(Cell.Types[i]),"gene_scaling.png",sep="")) 
    }
  }
}

# Gene Set Enrichment Analysis --------------------------------------------

#"Adaptive.Immunity"   
selected_genes <- Annotations %>% filter(Adaptive.Immunity == '+') %>% dplyr::select(c("Gene"),)

if(nrow(selected_genes)>=2){
  #filter genes by cell type
  #not all genes are within dataset, exclude those
  output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
  if(length(colnames(output))>1){
    output <- sapply(output, as.numeric)
    rownames(output) <- rownames(samples)
    
    #plot
    my_heatmap <- pheatmap::pheatmap(output,
                                     main="Adaptive.Immunity",
                                     annotation_row = samples,
                                     scale="column")
    #save heatmap
    save_pheatmap_png(my_heatmap, paste(as.character("Adaptive.Immunity"),"gene_scaling.png",sep="")) 
  }
}

#"Apoptosis"   
selected_genes <- Annotations %>% filter(Apoptosis == '+') %>% dplyr::select(c("Gene"),)

if(nrow(selected_genes)>=2){
  #filter genes by cell type
  #not all genes are within dataset, exclude those
  output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
  if(length(colnames(output))>1){
    output <- sapply(output, as.numeric)
    rownames(output) <- rownames(samples)
    
    #plot
    my_heatmap <- pheatmap::pheatmap(output,
                                     main="Apoptosis",
                                     annotation_row = samples,
                                     scale="column")
    #save heatmap
    save_pheatmap_png(my_heatmap, paste(as.character("Apoptosis"),"gene_scaling.png",sep="")) 
  }
}

#"Cell.Cycle"   
selected_genes <- Annotations %>% filter(Cell.Cycle == '+') %>% dplyr::select(c("Gene"),)

if(nrow(selected_genes)>=2){
  #filter genes by cell type
  #not all genes are within dataset, exclude those
  output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
  if(length(colnames(output))>1){
    output <- sapply(output, as.numeric)
    rownames(output) <- rownames(samples)
    
    #plot
    my_heatmap <- pheatmap::pheatmap(output,
                                     main="Cell.Cycle",
                                     annotation_row = samples,
                                     scale="column")
    #save heatmap
    save_pheatmap_png(my_heatmap, paste(as.character("Cell.Cycle"),"gene_scaling.png",sep="")) 
  }
}

#"Cellular.Stress"   
selected_genes <- Annotations %>% filter(Cellular.Stress == '+') %>% dplyr::select(c("Gene"),)

if(nrow(selected_genes)>=2){
  #filter genes by cell type
  #not all genes are within dataset, exclude those
  output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
  if(length(colnames(output))>1){
    output <- sapply(output, as.numeric)
    rownames(output) <- rownames(samples)
    
    #plot
    my_heatmap <- pheatmap::pheatmap(output,
                                     main="Cellular.Stress",
                                     annotation_row = samples,
                                     scale="column")
    #save heatmap
    save_pheatmap_png(my_heatmap, paste(as.character("Cellular.Stress"),"gene_scaling.png",sep="")) 
  }
}

#"Complement.System"   
selected_genes <- Annotations %>% filter(Complement.System == '+') %>% dplyr::select(c("Gene"),)

if(nrow(selected_genes)>=2){
  #filter genes by cell type
  #not all genes are within dataset, exclude those
  output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
  if(length(colnames(output))>1){
    output <- sapply(output, as.numeric)
    rownames(output) <- rownames(samples)
    
    #plot
    my_heatmap <- pheatmap::pheatmap(output,
                                     main="Complement.System",
                                     annotation_row = samples,
                                     scale="column")
    #save heatmap
    save_pheatmap_png(my_heatmap, paste(as.character("Complement.System"),"gene_scaling.png",sep="")) 
  }
}

#"Death.Receptor.Signaling"   
selected_genes <- Annotations %>% filter(Death.Receptor.Signaling == '+') %>% dplyr::select(c("Gene"),)

if(nrow(selected_genes)>=2){
  #filter genes by cell type
  #not all genes are within dataset, exclude those
  output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
  if(length(colnames(output))>1){
    output <- sapply(output, as.numeric)
    rownames(output) <- rownames(samples)
    
    #plot
    my_heatmap <- pheatmap::pheatmap(output,
                                     main="Death.Receptor.Signaling",
                                     annotation_row = samples,
                                     scale="column")
    #save heatmap
    save_pheatmap_png(my_heatmap, paste(as.character("Death.Receptor.Signaling"),"gene_scaling.png",sep="")) 
  }
}

#"Extracellular.matrix.organization"   
selected_genes <- Annotations %>% filter(Extracellular.matrix.organization == '+') %>% dplyr::select(c("Gene"),)

if(nrow(selected_genes)>=2){
  #filter genes by cell type
  #not all genes are within dataset, exclude those
  output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
  if(length(colnames(output))>1){
    output <- sapply(output, as.numeric)
    rownames(output) <- rownames(samples)
    
    #plot
    my_heatmap <- pheatmap::pheatmap(output,
                                     main="Extracellular.matrix.organization",
                                     annotation_row = samples,
                                     scale="column")
    #save heatmap
    save_pheatmap_png(my_heatmap, paste(as.character("Extracellular.matrix.organization"),"gene_scaling.png",sep="")) 
  }
}

#"Fc.Receptor.Signaling"   
selected_genes <- Annotations %>% filter(Fc.Receptor.Signaling == '+') %>% dplyr::select(c("Gene"),)

if(nrow(selected_genes)>=2){
  #filter genes by cell type
  #not all genes are within dataset, exclude those
  output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
  if(length(colnames(output))>1){
    output <- sapply(output, as.numeric)
    rownames(output) <- rownames(samples)
    
    #plot
    my_heatmap <- pheatmap::pheatmap(output,
                                     main="Fc.Receptor.Signaling",
                                     annotation_row = samples,
                                     scale="column")
    #save heatmap
    save_pheatmap_png(my_heatmap, paste(as.character("Fc.Receptor.Signaling"),"gene_scaling.png",sep="")) 
  }
}

#"Innate.Immunity"   
selected_genes <- Annotations %>% filter(Innate.Immunity == '+') %>% dplyr::select(c("Gene"),)

if(nrow(selected_genes)>=2){
  #filter genes by cell type
  #not all genes are within dataset, exclude those
  output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
  if(length(colnames(output))>1){
    output <- sapply(output, as.numeric)
    rownames(output) <- rownames(samples)
    
    #plot
    my_heatmap <- pheatmap::pheatmap(output,
                                     main="Innate.Immunity",
                                     annotation_row = samples,
                                     scale="column")
    #save heatmap
    save_pheatmap_png(my_heatmap, paste(as.character("Innate.Immunity"),"gene_scaling.png",sep="")) 
  }
}

#"Interferon.Signaling"   
selected_genes <- Annotations %>% filter(Interferon.Signaling == '+') %>% dplyr::select(c("Gene"),)

if(nrow(selected_genes)>=2){
  #filter genes by cell type
  #not all genes are within dataset, exclude those
  output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
  if(length(colnames(output))>1){
    output <- sapply(output, as.numeric)
    rownames(output) <- rownames(samples)
    
    #plot
    my_heatmap <- pheatmap::pheatmap(output,
                                     main="Interferon.Signaling",
                                     annotation_row = samples,
                                     scale="column")
    #save heatmap
    save_pheatmap_png(my_heatmap, paste(as.character("Interferon.Signaling"),"gene_scaling.png",sep="")) 
  }
}

#"Interleukin.Signaling"   
selected_genes <- Annotations %>% filter(Interleukin.Signaling == '+') %>% dplyr::select(c("Gene"),)

if(nrow(selected_genes)>=2){
  #filter genes by cell type
  #not all genes are within dataset, exclude those
  output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
  if(length(colnames(output))>1){
    output <- sapply(output, as.numeric)
    rownames(output) <- rownames(samples)
    
    #plot
    my_heatmap <- pheatmap::pheatmap(output,
                                     main="Interleukin.Signaling",
                                     annotation_row = samples,
                                     scale="column")
    #save heatmap
    save_pheatmap_png(my_heatmap, paste(as.character("Interleukin.Signaling"),"gene_scaling.png",sep="")) 
  }
}

#"NF.kB"   
selected_genes <- Annotations %>% filter(NF.kB == '+') %>% dplyr::select(c("Gene"),)

if(nrow(selected_genes)>=2){
  #filter genes by cell type
  #not all genes are within dataset, exclude those
  output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
  if(length(colnames(output))>1){
    output <- sapply(output, as.numeric)
    rownames(output) <- rownames(samples)
    
    #plot
    my_heatmap <- pheatmap::pheatmap(output,
                                     main="NF.kB",
                                     annotation_row = samples,
                                     scale="column")
    #save heatmap
    save_pheatmap_png(my_heatmap, paste(as.character("NF.kB"),"gene_scaling.png",sep="")) 
  }
}

#"Metabolism"   
selected_genes <- Annotations %>% filter(Metabolism == '+') %>% dplyr::select(c("Gene"),)

if(nrow(selected_genes)>=2){
  #filter genes by cell type
  #not all genes are within dataset, exclude those
  output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
  if(length(colnames(output))>1){
    output <- sapply(output, as.numeric)
    rownames(output) <- rownames(samples)
    
    #plot
    my_heatmap <- pheatmap::pheatmap(output,
                                     main="Metabolism",
                                     annotation_row = samples,
                                     scale="column")
    #save heatmap
    save_pheatmap_png(my_heatmap, paste(as.character("Metabolism"),"gene_scaling.png",sep="")) 
  }
}

#"TLR.Signaling"   
selected_genes <- Annotations %>% filter(TLR.Signaling == '+') %>% dplyr::select(c("Gene"),)

if(nrow(selected_genes)>=2){
  #filter genes by cell type
  #not all genes are within dataset, exclude those
  output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
  if(length(colnames(output))>1){
    output <- sapply(output, as.numeric)
    rownames(output) <- rownames(samples)
    
    #plot
    my_heatmap <- pheatmap::pheatmap(output,
                                     main="TLR.Signaling",
                                     annotation_row = samples,
                                     scale="column")
    #save heatmap
    save_pheatmap_png(my_heatmap, paste(as.character("TLR.Signaling"),"gene_scaling.png",sep="")) 
  }
}

#"VEGF.Signaling"   
selected_genes <- Annotations %>% filter(VEGF.Signaling == '+') %>% dplyr::select(c("Gene"),)

if(nrow(selected_genes)>=2){
  #filter genes by cell type
  #not all genes are within dataset, exclude those
  output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
  if(length(colnames(output))>1){
    output <- sapply(output, as.numeric)
    rownames(output) <- rownames(samples)
    
    #plot
    my_heatmap <- pheatmap::pheatmap(output,
                                     main="VEGF.Signaling",
                                     annotation_row = samples,
                                     scale="column")
    #save heatmap
    save_pheatmap_png(my_heatmap, paste(as.character("VEGF.Signaling"),"gene_scaling.png",sep="")) 
  }}
    
#"Wnt.Signaling"   
selected_genes <- Annotations %>% filter(Wnt.Signaling == '+') %>% dplyr::select(c("Gene"),)
    
if(nrow(selected_genes)>=2){
  #filter genes by cell type
  #not all genes are within dataset, exclude those
  output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
  if(length(colnames(output))>1){
    output <- sapply(output, as.numeric)
    rownames(output) <- rownames(samples)
    
    #plot
    my_heatmap <- pheatmap::pheatmap(output,
                                     main="Wnt.Signaling",
                                     annotation_row = samples,
                                     scale="column")
    #save heatmap
    save_pheatmap_png(my_heatmap, paste(as.character("Wnt.Signaling"),"gene_scaling.png",sep="")) 
  }}





# Misc --------------------------------------------------------------------

# Immune Cell Type Composition?

What are the expression heatmaps across genes that are suppose to indicate cell type calls

### Heatmaps of Each Immune Cell Type

```{r,echo=FALSE}
#rotate expression data so genes are columns and clean it up
data <- as_tibble(cbind(nms = names(Expression_Data), t(Expression_Data)))
colnames(data) <-data[2,]
data <- data[-1, ]
data <- data[-1, ]
data <- data[-1, ]



#meta information orientated to data
samples <- data.frame(data$Gene_Name)
colnames(samples) <- "Sample.Name"
samples <- samples %>% left_join(Metadata)
rownames(samples) <- samples$Sample.Name
samples <- samples %>% dplyr::select(-c("Sample.Name","File.Name","Cartridge.ID","Lane.Number","Import.Date")) #remove unncessary colnames

#function for saving
save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


#Cell Types Available
Cell.Types <- unique(Annotations$Cell.Type)

for(i in 2:length(Cell.Types)){
  #Select genes available for Cell Type
  selected_genes <- Annotations %>% filter(Cell.Type == Cell.Types[i])
  
  if(nrow(selected_genes)>=2){
    #filter genes by cell type
    #not all genes are within dataset, exclude those
    output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
    if(length(colnames(output))>1){
      output <- sapply(output, as.numeric)
      rownames(output) <- data$Gene_Name
      
      #plot
      my_heatmap <- pheatmap::pheatmap(output,main=Cell.Types[i],annotation_row = samples,scale="column")
      #save heatmap
      save_pheatmap_png(my_heatmap, paste(as.character(Cell.Types[i]),"gene_scaling.png",sep="")) 
    }
  }
}

```

### Excluding Samples with Treatment 

```{r,echo=FALSE}

#rotate expression data so genes are columns and clean it up
data <- as_tibble(cbind(nms = names(Expression_Data), t(Expression_Data)))
colnames(data) <-data[2,]
data <- data[-1, ]
data <- data[-1, ]
data <- data[-1, ]

#exclude samples without treatment
data <- data %>% filter(Gene_Name %in% ((Metadata %>% filter(Treatment == 'None'))$Sample.Name %>% as.character()))

#meta information orientated to data
samples <- data.frame(data$Gene_Name)
colnames(samples) <- "Sample.Name"
samples <- samples %>% left_join(Metadata)
rownames(samples) <- samples$Sample.Name
samples <- samples %>% dplyr::select(-c("Sample.Name","File.Name","Cartridge.ID","Lane.Number","Import.Date")) #remove unncessary colnames

#function for saving
save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


#Cell Types Available
Cell.Types <- unique(Annotations$Cell.Type)

for(i in 2:length(Cell.Types)){
  #Select genes available for Cell Type
  selected_genes <- Annotations %>% filter(Cell.Type == Cell.Types[i])
  
  if(nrow(selected_genes)>=2){
    #filter genes by cell type
    #not all genes are within dataset, exclude those
    output <- data %>% dplyr::select(any_of(selected_genes$Gene)) 
    if(length(colnames(output))>1){
      output <- sapply(output, as.numeric)
      rownames(output) <- data$Gene_Name
      
      #plot
      my_heatmap <- pheatmap::pheatmap(output,main=Cell.Types[i],annotation_row = samples,scale="column")
      #save heatmap
      save_pheatmap_png(my_heatmap, paste(as.character(Cell.Types[i]),"gene_scaling.png",sep="")) 
    }
  }
}

```
