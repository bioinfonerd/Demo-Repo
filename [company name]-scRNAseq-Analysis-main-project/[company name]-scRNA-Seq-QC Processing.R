# Purpose: Complete scRNA-Seq Analysis on [company name] scRNA-Seq Data


# Libraries ---------------------------------------------------------------

# required ubuntu server files
# To install
# sudo apt install r-base-dev libfontconfig1-dev pkgconf libharfbuzz-dev
# libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
# cmake gsl-bin libgsl-dev cairo-5c libcairo-5c-dev libcairo-5c0
# libxt-dev libcairo2-dev libhdf5-dev

#Use renv to organize R library as local instant copy
#Setup R project, enable renv for environment organization
#alternatively as renv sometimes does not load correctly, just update server libraries

##############################################################
# To efficiently handle installation and loading of Libraries#
##############################################################
# Package names (rcran packages)
packages <- c("ggplot2", "tidyr", "dplyr", "Seurat", "tidyverse", "knitr",
              "Matrix", "reshape2", "ggExtra", "patchwork","scales",
              "cowplot", "circlize", "ggridges",
              "ggpmisc", "glmnet", "foreach","ggrepel", "EnvStats", "corrr",
              "igraph", "ggraph","ggplot.multistats", "harmony",
              "ggpubr", "survminer", "data.table",
              "kableExtra", "SCINA", "SAVER", "future","sctransform")

# Install packages not yet installed (rcran packages)
# install.packages('devtools')
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], dependencies = TRUE)
}

# Package names (Bioconductor packages)
bioconductor_packages <- c("SingleCellExperiment", "SingleR", "ComplexHeatmap", 
              "enrichplot","scater", "AnnotationHub","clusterProfiler",
              "scRepertoire","glmGamPoi")

# Install packages not yet installed 
# install.packages("BiocManager")
installed_bioconductor_packages <- bioconductor_packages %in% rownames(installed.packages())
if (any(installed_bioconductor_packages == FALSE)) {
  BiocManager::install(bioconductor_packages[!installed_bioconductor_packages])
}

# Github Package Installs
github_packages <- c("SeuratWrappers","SeuratDisk","SeuratData","ProjecTILs")

#remotes::install_github('satijalab/seurat-wrappers')
#remotes::install_github("mojaveazure/seurat-disk")
#devtools::install_github('satijalab/seurat-data')
#remotes::install_github("carmonalab/STACAS")
#remotes::install_github("carmonalab/ProjecTILs")

#no longer available
# Rmagic

# devtools::install_github("ncborcherding/scRepertoire@dev")

# Library loading
invisible(lapply(packages, library, character.only = TRUE))
invisible(lapply(bioconductor_packages, library, character.only = TRUE))
invisible(lapply(github_packages, library, character.only = TRUE))

# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# Global Variables --------------------------------------------------------

#setup experiment conditions
experiment_name = "[company name] Example"
dataset_loc <- "/home/njohnson/data/aligned_data"
all_ids <- c("SE2_Drug_1","SE2_blemocyin","SE3_Drug_A",
         "SE3_bleo","_SE2_Drug_3A","SE2_Drug_2","SE2_control",
         "SE3_Drug_B","SE3_control")

se2_ids <- c("SE2_Drug_1","SE2_blemocyin","_SE2_Drug_3A","SE2_Drug_2","SE2_control")

se3_ids <- c("SE3_Drug_A","SE3_bleo","SE3_Drug_B","SE3_control")



# Load QC Summary Data & Generate Summary Table---------------------------------

#Metrics
d10x.metrics <- lapply(ids, function(i){
  # remove _Counts is if names don't include them
  metrics <- read.csv(file.path(dataset_loc,paste0(i,"/outs"),"metrics_summary.csv"), colClasses = "character")
})
experiment.metrics <- do.call("rbind", d10x.metrics)
rownames(experiment.metrics) <- ids

sequencing_metrics <- data.frame(t(experiment.metrics[,c(4:16,1,17,2,3,18,19)]))
row.names(sequencing_metrics) <- gsub("\\."," ", rownames(sequencing_metrics))

#create table of data
sequencing_metrics %>%
  kable(caption = 'Cell Ranger Results') %>%
  pack_rows("Sequencing Characteristics", 1, 6, label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Mapping Characteristics", 7, 13, label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Cell Characteristics", 14, 19, label_row_css = "background-color: #666; color: #fff;") %>%
  kable_styling("striped")


# Create Seurat Object ----------------------------------------------------

#load in each experiment data and all data to create separate Seurat objects 
# (uses alot of RAM, best to work with each object and save-current instance is large enough so won't separate)

#load in data
d10x.all.data <- lapply(all_ids, function(i){
  d10x <- Read10X_h5(file.path(dataset_loc,paste0(i,"/outs"),"raw_feature_bc_matrix.h5"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})
names(d10x.all.data) <- all_ids

d10x.se2.data <- lapply(se2_ids, function(i){
  d10x <- Read10X_h5(file.path(dataset_loc,paste0(i,"/outs"),"raw_feature_bc_matrix.h5"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})
names(d10x.se2.data) <- se2_ids

d10x.se3.data <- lapply(se3_ids, function(i){
  d10x <- Read10X_h5(file.path(dataset_loc,paste0(i,"/outs"),"raw_feature_bc_matrix.h5"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})
names(d10x.se3.data) <- se3_ids

#combine dataset
experiment.all.data <- do.call("cbind", d10x.all.data)
experiment.se2.data <- do.call("cbind", d10x.se2.data)
experiment.se3.data <- do.call("cbind", d10x.se3.data)

#create Seurat Object
#filter criteria: 
#remove genes that do not occur in a minimum of 0 cells
#remove cells that don’t have a minimum of 200 features

experiment.all.aggregate.filtered <- CreateSeuratObject(
  experiment.all.data,
  project = "All Samples",
  min.cells = 1,
  min.features = 300,
  names.field = 2,
  names.delim = "\\-")

experiment.se3.aggregate.filtered <- CreateSeuratObject(
  experiment.se3.data,
  project = "SE3 Experiment",
  min.cells = 1,
  min.features = 300,
  names.field = 2,
  names.delim = "\\-")

experiment.se2.aggregate.filtered <- CreateSeuratObject(
  experiment.se2.data,
  project = "SE2 Experiment",
  min.cells = 1,
  min.features = 300,
  names.field = 2,
  names.delim = "\\-")

#MIT genes
experiment.all.aggregate.filtered[["percent.mt"]] <- PercentageFeatureSet(experiment.all.aggregate.filtered, pattern = "^MT-")
experiment.se3.aggregate.filtered[["percent.mt"]] <- PercentageFeatureSet(experiment.se3.aggregate.filtered, pattern = "^MT-")
experiment.se2.aggregate.filtered[["percent.mt"]] <- PercentageFeatureSet(experiment.se2.aggregate.filtered, pattern = "^MT-")

#remove unused objects
rm(d10x.all.data)
rm(d10x.se3.data)
rm(d10x.se2.data)

rm(experiment.all.data)
rm(experiment.se2.data)
rm(experiment.se3.data)

#save object
saveRDS(experiment.all.aggregate.filtered, file = "experiment.all.aggregate.filtered-fresh.rds")
saveRDS(experiment.se3.aggregate.filtered, file = "experiment.se3.aggregate.filtered-fresh.rds")
saveRDS(experiment.se2.aggregate.filtered, file = "experiment.se2.aggregate.filtered-fresh.rds")

# Run QC Plots -----------------------------------------------------------------

#(must be done before seurat object is updated)
# will overwrite sample id names

#All
png(file = "All QC Metrics - percent mt.png", width = 600, height = 600) 
VlnPlot(experiment.all.aggregate.filtered,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3,
        pt.size = 0.2,
        split.by= "orig.ident")
dev.off()

png(file = "All QC Metrics - wo percent mt.png", width = 600, height = 600) 
VlnPlot(experiment.all.aggregate.filtered, 
        features = c("nFeature_RNA", "nCount_RNA"),
        ncol = 2,
        pt.size = 0.2,
        split.by= "orig.ident")
dev.off()

#SE2
png(file = "SE2 QC Metrics - percent mt.png", width = 600, height = 600) 
VlnPlot(experiment.se2.aggregate.filtered,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3,
        pt.size = 0.2,
        split.by= "orig.ident")
dev.off()

png(file = "SE2 QC Metrics - wo percent mt.png", width = 600, height = 600) 
VlnPlot(experiment.se2.aggregate.filtered, 
        features = c("nFeature_RNA", "nCount_RNA"),
        ncol = 2,
        pt.size = 0.2,
        split.by= "orig.ident")
dev.off()


#SE3
png(file = "SE3 QC Metrics - percent mt.png", width = 600, height = 600) 
VlnPlot(experiment.se3.aggregate.filtered,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3,
        pt.size = 0.2,
        split.by= "orig.ident")
dev.off()

png(file = "SE3 QC Metrics - wo percent mt.png", width = 600, height = 600) 
VlnPlot(experiment.se3.aggregate.filtered, 
        features = c("nFeature_RNA", "nCount_RNA"),
        ncol = 2,
        pt.size = 0.2,
        split.by= "orig.ident")
dev.off()




# Normalize Data ---------------------------------------------------------------

#ALL data
experiment.all.aggregate.filtered <- SCTransform(experiment.all.aggregate.filtered,
                                             method = "glmGamPoi",
                                             verbose = TRUE)

#SE2 data
experiment.se2.aggregate.filtered <- SCTransform(experiment.se2.aggregate.filtered,
                                                 method = "glmGamPoi",
                                                 verbose = TRUE)

#SE3 data
experiment.se3.aggregate.filtered <- SCTransform(experiment.se3.aggregate.filtered,
                                                 method = "glmGamPoi",
                                                 verbose = TRUE)

saveRDS(experiment.all.aggregate.filtered, file = "experiment.all.aggregate.filtered-normalized.rds")
saveRDS(experiment.se2.aggregate.filtered, file = "experiment.se2.aggregate.filtered-normalized.rds")
saveRDS(experiment.se3.aggregate.filtered, file = "experiment.se3.aggregate.filtered-normalized.rds")

# Identify Highest Variable Genes -----------------------------------------

#Find Variable Features
experiment.all.aggregate.filtered <- FindVariableFeatures(experiment.all.aggregate.filtered, selection.method = "vst", nfeatures = 2000)
experiment.se2.aggregate.filtered <- FindVariableFeatures(experiment.se2.aggregate.filtered, selection.method = "vst", nfeatures = 2000)
experiment.se3.aggregate.filtered <- FindVariableFeatures(experiment.se3.aggregate.filtered, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
all.top10 <- head(VariableFeatures(experiment.all.aggregate.filtered), 10)
se2.top10 <- head(VariableFeatures(experiment.se2.aggregate.filtered), 10)
se3.top10 <- head(VariableFeatures(experiment.se3.aggregate.filtered), 10)

# plot variable features with labels

png(file = "Variable Genes - ALL.png", width = 600, height = 600) 
plot1 <- VariableFeaturePlot(experiment.all.aggregate.filtered)
LabelPoints(plot = plot1, points = all.top10, repel = TRUE)
dev.off()

png(file = "Variable Genes - SE2.png", width = 600, height = 600) 
plot1 <- VariableFeaturePlot(experiment.se2.aggregate.filtered)
LabelPoints(plot = plot1, points = se2.top10, repel = TRUE)
dev.off()

png(file = "Variable Genes - SE3.png", width = 600, height = 600) 
plot1 <- VariableFeaturePlot(experiment.se3.aggregate.filtered)
LabelPoints(plot = plot1, points = se3.top10, repel = TRUE)
dev.off()

# Dimension Reduction & Clustering ---------------------------------------------

# All Samples
experiment.all.aggregate.filtered <- RunPCA(experiment.all.aggregate.filtered, verbose = TRUE)
experiment.all.aggregate.filtered <- RunUMAP(experiment.all.aggregate.filtered, dims = 1:30, verbose = TRUE)
experiment.all.aggregate.filtered <- FindNeighbors(experiment.all.aggregate.filtered, dims = 1:30, verbose = TRUE)
experiment.all.aggregate.filtered <- FindClusters(experiment.all.aggregate.filtered, verbose = TRUE)

# SE2
experiment.se2.aggregate.filtered <- RunPCA(experiment.se2.aggregate.filtered, verbose = TRUE)
experiment.se2.aggregate.filtered <- RunUMAP(experiment.se2.aggregate.filtered, dims = 1:30, verbose = TRUE)
experiment.se2.aggregate.filtered <- FindNeighbors(experiment.se2.aggregate.filtered, dims = 1:30, verbose = TRUE)
experiment.se2.aggregate.filtered <- FindClusters(experiment.se2.aggregate.filtered, verbose = TRUE)

# SE3 
experiment.se3.aggregate.filtered <- RunPCA(experiment.se3.aggregate.filtered, verbose = TRUE)
experiment.se3.aggregate.filtered <- RunUMAP(experiment.se3.aggregate.filtered, dims = 1:30, verbose = TRUE)
experiment.se3.aggregate.filtered <- FindNeighbors(experiment.se3.aggregate.filtered, dims = 1:30, verbose = TRUE)
experiment.se3.aggregate.filtered <- FindClusters(experiment.se3.aggregate.filtered, verbose = TRUE)

#Plot UMAP + Clusters

#All Samples
png(file = "Sample-Umap-Cluster-ALL samples-Combined.png", width = 600, height = 600)
DimPlot(experiment.all.aggregate.filtered,
        label = TRUE,
        raster = FALSE) + NoLegend()
dev.off()

png(file = "Sample-Umap-Cluster-ALL samples-By Sample.png", width = 600, height = 600)
#Plot by sample type UMAP + clustering
DimPlot(experiment.all.aggregate.filtered,
        label = TRUE,
        split.by= "orig.ident",
        raster = FALSE) + NoLegend()
dev.off()

#SE2 
png(file = "Sample-Umap-Cluster-SE2 samples-Combined.png", width = 600, height = 600)
DimPlot(experiment.se2.aggregate.filtered,
        label = TRUE,
        raster = FALSE) + NoLegend()
dev.off()

png(file = "Sample-Umap-Cluster-SE2 samples-By Sample.png", width = 600, height = 600)
#Plot by sample type UMAP + clustering
DimPlot(experiment.se2.aggregate.filtered,
        label = TRUE,
        split.by= "orig.ident",
        raster = FALSE) + NoLegend()
dev.off()

#SE3
png(file = "Sample-Umap-Cluster-SE3 samples-Combined.png", width = 600, height = 600)
DimPlot(experiment.se3.aggregate.filtered,
        label = TRUE,
        raster = FALSE) + NoLegend()
dev.off()

png(file = "Sample-Umap-Cluster-SE3 samples-By Sample.png", width = 600, height = 600)
#Plot by sample type UMAP + clustering
DimPlot(experiment.se3.aggregate.filtered,
        label = TRUE,
        split.by= "orig.ident",
        raster = FALSE) + NoLegend()
dev.off()

#save data
saveRDS(experiment.all.aggregate.filtered, file = "experiment.all.aggregate.filtered-normalized-dim-cluster.rds")
saveRDS(experiment.se2.aggregate.filtered, file = "experiment.se2.aggregate.filtered-normalized-dim-cluster.rds")
saveRDS(experiment.se3.aggregate.filtered, file = "experiment.se3.aggregate.filtered-normalized-dim-cluster.rds")

# Cluster Biomarkers -----------------------------------------------------

# find markers for every cluster compared to all remaining cells, report only the positive ones

experiment.all.aggregate.filtered.cluster.markers <- FindAllMarkers(experiment.all.aggregate.filtered,
                                                                only.pos = TRUE,
                                                                min.pct = 0.25,
                                                                logfc.threshold = 0.25)

experiment.se2.aggregate.filtered.cluster.markers <- FindAllMarkers(experiment.se2.aggregate.filtered,
                                                                only.pos = TRUE,
                                                                min.pct = 0.25,
                                                                logfc.threshold = 0.25)

experiment.se3.aggregate.filtered.cluster.markers <- FindAllMarkers(experiment.se3.aggregate.filtered,
                                                                only.pos = TRUE,
                                                                min.pct = 0.25,
                                                                logfc.threshold = 0.25)

# save results
saveRDS(experiment.all.aggregate.filtered.cluster.markers, file = "experiment.all.aggregate.filtered.cluster.markers.rds")
saveRDS(experiment.se2.aggregate.filtered.cluster.markers, file = "experiment.se2.aggregate.filtered.cluster.markers.rds")
saveRDS(experiment.se3.aggregate.filtered.cluster.markers, file = "experiment.se3.aggregate.filtered.cluster.markers.rds")



# Biomarker Cluster Analysis ----------------------------------------------

experiment.aggregate.filtered.cluster.markers.base.se2_control <- FindAllMarkers(experiment.aggregate.filtered,
                                                                                 only.pos = TRUE,
                                                                                 min.pct = 0.25,
                                                                                 logfc.threshold = 0.25)

experiment.aggregate.filtered.cluster.markers.base.se3_control %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

experiment.aggregate.filtered.cluster.markers.base.se2_control %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)


experiment.aggregate.filtered.cluster.markers.base.se3_control %>%
  group_by(cluster) %>%
  arrange(avg_log2FC) %>%
  filter(gene == "Col1a1")

experiment.aggregate.filtered.cluster.markers.base.se2_control %>%
  group_by(cluster) %>%
  arrange(avg_log2FC) %>%
  filter(gene == "Col1a1")

experiment.aggregate.filtered.cluster.markers %>%
  filter(cluster == 18) %>%
  arrange(p_val_adj)

load("experiment.aggregate.filtered.cluster.markers.rds")


# Plotting Genes of Interest for each sample as violin plots--------------------

# Visualize canonical marker genes as violin plots.

# Profibrotic Genes
png(file = "Profibrotic Genes-across dataset.png", width = 600, height = 600) 
VlnPlot(experiment.aggregate.filtered,
        features = c("Col1a1", "Acta2", "Timp1", "Serpine1", "Fn1"),
        pt.size=0,
        ncol = 3) 
dev.off()

#Senescence related genes
png(file = "Senescence related genes-across dataset.png") 
VlnPlot(experiment.aggregate.filtered,
        features = c("Cdkn1a", "Cdkn1b", "Cdkn1a",
                     "Cdkn2b", "Glb1","Tp53","Serpine1"),
        pt.size=0,
        ncol = 3)
dev.off()

#	SASP genes (secretory associated senescence phenotype)
png(file = "SASP genes-across dataset.png") 
VlnPlot(experiment.aggregate.filtered,
        features = c("Il6", "Cxcl8", "Ccl2",
                     "Vegf", "Fgf2","Tgfb1",
                     "Egf","Cxcl1","Col1a1"),
        pt.size=0,
        ncol = 3)
dev.off()

#	Ephrin receptor and ligands genes
# Ephrin Ligands
png(file = "Ephrin Ligand genes-across dataset.png") 
VlnPlot(experiment.aggregate.filtered,
        features = c("Efna1", "Efna2", "Efna3",
                     "Efna4", "Efna5","Efnb1",
                     "Efnb2","Efnb3"),
        pt.size=0,
        ncol = 3)
dev.off()

# Ephrin receptor
png(file = "Ephrin Receptor genes-across dataset.png") 
VlnPlot(experiment.aggregate.filtered,
        features = c("Epha1", "Epha2", "Epha3",
                     "Epha4", "Epha5","Epha6",
                     "Epha7","Epha8","Epha10",
                     "Ephb1","Ephb2","Ephb3",
                     "Ephb4","Ephb6"),
        pt.size=0,
        ncol = 3)
dev.off()

# BCL2 family of proteins 
png(file = "BCL2 genes-across dataset.png")
VlnPlot(experiment.aggregate.filtered,
        features = c("Bax", "Bak1", "Bok",
                     "Bcl2l11", "Bid","Bcl2",
                     "Bcl2l2","Bcl2l1","Bcl2a1",
                     "Mcl1","Bcl2l10","Bbc3",
                     "Bmf","Bad","Pmaip1","Hrk",
                     "Bnip3","Bik"),
        pt.size=0,
        ncol = 3)
dev.off()

#Significant genes [company name] Ephrin treatment:
# Condition 474
png(file = "Condition 474-across dataset.png")
VlnPlot(experiment.aggregate.filtered,
        features = c("Gm5869", "Nr1d1", "Abtb1",
                     "Dbp", "Prdm8","Rpl3-ps2",
                     "1110038F14Rik","Nfil3","H2-Q6",
                     "Sppl3","Bhlhe41","Nr4a3",
                     "Gm22918","Ccdc85b","Pip4p1","Sema6c",
                     "Gm8229","Gata3","Mcpt4","Tmem252"),
        pt.size=0,
        ncol = 3)
dev.off()

# Condition 124 (B11)
png(file = "Condition 124-across dataset.png")
VlnPlot(experiment.aggregate.filtered,
        features = c("Rbm25", "Gm11263", "Msl3l2",
                     "Rnf187", "Stub1","Nr1d1",
                     "Nub1","Gm5869","Cops8",
                     "Oaz1-ps","Tmem41a","Mtch1",
                     "Rpl3-ps2","Thap11","Kif20b","Chpf",
                     "Igkv3-2","Ppp2cb","Emc10","Zfp574"),
        pt.size=0,
        ncol = 3)
dev.off()

#	Condition 394
png(file = "Condition 394-across dataset.png")
VlnPlot(experiment.aggregate.filtered,
        features = c("Nr1d1", "Tmem252", "Camkk1",
                     "Gm11263", "Plekhf1","Id2",
                     "Nr4a3","Msl3l2","Bhlhe41",
                     "Adm","Nfil3","Cebpd",
                     "Kctd11","Socs3","Tmem8","Slc38a2",
                     "1200007C13Rik","Dbp","Myc","4930516B21Rik"),
        pt.size=0,
        ncol = 3)
dev.off()

#	Condition 470
png(file = "Condition 470-across dataset.png")
VlnPlot(experiment.aggregate.filtered,
        features = c("Arntl", "Dbp", "Irf2bp2",
                     "Oaz1-ps", "Hspa5","Bhlhe41",
                     "Npas2","Tmem252","Pip4p1",
                     "Hspa8","Adm","Nfil3",
                     "Stub1","Ywhaz","Rrbp1","Sec31a",
                     "Tomm20","Rasl11a","Id2","Rab17"),
        pt.size=0,
        ncol = 3)
dev.off()

# Visualization on Umap ---------------------------------------------------

# Profibrotic Genes
png(file = "Profibrotic Genes-UMAP.png", width = 600, height = 600) 
FeaturePlot(experiment.aggregate.filtered,
        features = c("Col1a1", "Acta2", "Timp1", "Serpine1", "Fn1"),
        pt.size = 0.2,
        ncol = 3) 
dev.off()

#Senescence related genes
png(file = "Senescence related genes-UMAP.png") 
FeaturePlot(experiment.aggregate.filtered,
        features = c("Cdkn1a", "Cdkn1b", "Cdkn1a",
                     "Cdkn2b", "Glb1","Tp53","Serpine1"),
        pt.size=0.2,
        ncol = 3)
dev.off()

#	SASP genes (secretory associated senescence phenotype)
png(file = "SASP genes-UMAP.png") 
FeaturePlot(experiment.aggregate.filtered,
        features = c("Il6", "Cxcl8", "Ccl2",
                     "Vegf", "Fgf2","Tgfb1",
                     "Egf","Cxcl1","Col1a1"),
        pt.size=0.2,
        ncol = 3)
dev.off()

#	Ephrin receptor and ligands genes
# Ephrin Ligands
png(file = "Ephrin Ligand genes-UMAP.png") 
FeaturePlot(experiment.aggregate.filtered,
        features = c("Efna1", "Efna2", "Efna3",
                     "Efna4", "Efna5","Efnb1",
                     "Efnb2","Efnb3"),
        pt.size=0.2,
        ncol = 3)
dev.off()

# Ephrin receptor
png(file = "Ephrin Receptor genes-UMAP.png") 
FeaturePlot(experiment.aggregate.filtered,
        features = c("Epha1", "Epha2", "Epha3",
                     "Epha4", "Epha5","Epha6",
                     "Epha7","Epha8","Epha10",
                     "Ephb1","Ephb2","Ephb3",
                     "Ephb4","Ephb6"),
        pt.size=0.2,
        ncol = 3)
dev.off()

# BCL2 family of proteins 
png(file = "BCL2 genes-UMAP.png")
FeaturePlot(experiment.aggregate.filtered,
        features = c("Bax", "Bak1", "Bok",
                     "Bcl2l11", "Bid","Bcl2",
                     "Bcl2l2","Bcl2l1","Bcl2a1",
                     "Mcl1","Bcl2l10","Bbc3",
                     "Bmf","Bad","Pmaip1","Hrk",
                     "Bnip3","Bik"),
        pt.size=0.2,
        ncol = 3)
dev.off()

#Significant genes [company name] Ephrin treatment:
# Condition 474
png(file = "Condition 474-UMAP.png")
FeaturePlot(experiment.aggregate.filtered,
        features = c("Gm5869", "Nr1d1", "Abtb1",
                     "Dbp", "Prdm8","Rpl3-ps2",
                     "1110038F14Rik","Nfil3","H2-Q6",
                     "Sppl3","Bhlhe41","Nr4a3",
                     "Gm22918","Ccdc85b","Pip4p1","Sema6c",
                     "Gm8229","Gata3","Mcpt4","Tmem252"),
        pt.size=0.2,
        ncol = 3)
dev.off()

# Condition 124 (B11)
png(file = "Condition 124-UMAP.png")
FeaturePlot(experiment.aggregate.filtered,
        features = c("Rbm25", "Gm11263", "Msl3l2",
                     "Rnf187", "Stub1","Nr1d1",
                     "Nub1","Gm5869","Cops8",
                     "Oaz1-ps","Tmem41a","Mtch1",
                     "Rpl3-ps2","Thap11","Kif20b","Chpf",
                     "Igkv3-2","Ppp2cb","Emc10","Zfp574"),
        pt.size=0.2,
        ncol = 3)
dev.off()

#	Condition 394
png(file = "Condition 394-UMAP.png")
FeaturePlot(experiment.aggregate.filtered,
        features = c("Nr1d1", "Tmem252", "Camkk1",
                     "Gm11263", "Plekhf1","Id2",
                     "Nr4a3","Msl3l2","Bhlhe41",
                     "Adm","Nfil3","Cebpd",
                     "Kctd11","Socs3","Tmem8","Slc38a2",
                     "1200007C13Rik","Dbp","Myc","4930516B21Rik"),
        pt.size=0.2,
        ncol = 3)
dev.off()

#	Condition 470
png(file = "Condition 470-UMAP.png")
FeaturePlot(experiment.aggregate.filtered,
        features = c("Arntl", "Dbp", "Irf2bp2",
                     "Oaz1-ps", "Hspa5","Bhlhe41",
                     "Npas2","Tmem252","Pip4p1",
                     "Hspa8","Adm","Nfil3",
                     "Stub1","Ywhaz","Rrbp1","Sec31a",
                     "Tomm20","Rasl11a","Id2","Rab17"),
        pt.size=0.2,
        ncol = 3)
dev.off()





# Gene Focused ------------------------------------------------------------

# Sftpc
png(file = "Sftpc-across ALL dataset.png", width = 600, height = 600) 
VlnPlot(`experiment.all.aggregate.filtered-normalized-dim-cluster`,
        features = c("Sftpc"),
        pt.size=0,
        ncol = 1) 
dev.off()

png(file = "Sftpc-across ALL UMAP.png", width = 600, height = 600) 
FeaturePlot(`experiment.all.aggregate.filtered-normalized-dim-cluster`,
            features = c("Sftpc"),
            pt.size = 0.2,
            ncol = 1) 
dev.off()


png(file = "Sftpc-across SE2 dataset.png", width = 600, height = 600) 
VlnPlot(`experiment.se2.aggregate.filtered-normalized-dim-cluster`,
        features = c("Sftpc"),
        pt.size=0,
        ncol = 1) 
dev.off()

png(file = "Sftpc-across SE2 UMAP.png", width = 600, height = 600) 
FeaturePlot(`experiment.se2.aggregate.filtered-normalized-dim-cluster`,
            features = c("Sftpc"),
            pt.size = 0.2,
            ncol = 1) 
dev.off()


png(file = "Sftpc-across SE3 dataset.png", width = 600, height = 600) 
VlnPlot(`experiment.se3.aggregate.filtered-normalized-dim-cluster`,
        features = c("Sftpc"),
        pt.size=0,
        ncol = 1) 
dev.off()

png(file = "Sftpc-across SE3 UMAP.png", width = 600, height = 600) 
FeaturePlot(`experiment.se3.aggregate.filtered-normalized-dim-cluster`,
            features = c("Sftpc"),
            pt.size = 0.2,
            ncol = 1) 
dev.off()






# Lyz2
png(file = "Lyz2-across ALL dataset.png", width = 600, height = 600) 
VlnPlot(`experiment.all.aggregate.filtered-normalized-dim-cluster`,
        features = c("Lyz2"),
        pt.size=0,
        ncol = 1) 
dev.off()

png(file = "Lyz2-across ALL UMAP.png", width = 600, height = 600) 
FeaturePlot(`experiment.all.aggregate.filtered-normalized-dim-cluster`,
            features = c("Lyz2"),
            pt.size = 0.2,
            ncol = 1) 
dev.off()


png(file = "Lyz2-across SE2 dataset.png", width = 600, height = 600) 
VlnPlot(`experiment.se2.aggregate.filtered-normalized-dim-cluster`,
        features = c("Lyz2"),
        pt.size=0,
        ncol = 1) 
dev.off()

png(file = "Lyz2-across SE2 UMAP.png", width = 600, height = 600) 
FeaturePlot(`experiment.se2.aggregate.filtered-normalized-dim-cluster`,
            features = c("Lyz2"),
            pt.size = 0.2,
            ncol = 1) 
dev.off()


png(file = "Lyz2-across SE3 dataset.png", width = 600, height = 600) 
VlnPlot(`experiment.se3.aggregate.filtered-normalized-dim-cluster`,
        features = c("Lyz2"),
        pt.size=0,
        ncol = 1) 
dev.off()

png(file = "Lyz2-across SE3 UMAP.png", width = 600, height = 600) 
FeaturePlot(`experiment.se3.aggregate.filtered-normalized-dim-cluster`,
            features = c("Lyz2"),
            pt.size = 0.2,
            ncol = 1) 
dev.off()










# Cell Type Annotation (scAnnotate) ---------------------------------------
# https://academic.oup.com/bioinformaticsadvances/article/3/1/vbad030/7076619


# Cell Type Annotation (SingleR) ------------------------------------------
# https://bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html
# https://bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.R
# https://bioconductor.org/packages/release/bioc/manuals/SingleR/man/SingleR.pdf

# Cell Type Annotation (ScType) ----------------------------------
#https://www.nature.com/articles/s41467-022-28803-w, https://github.com/IanevskiAleksandr/sc-type/

#SE2
experiment.se2.aggregate.filtered.normalized.dim.cluster <- readRDS("~/data/R-scRNA-Seq-Processing/R Objects/experiment.se2.aggregate.filtered-normalized-dim-cluster.rds")

# Cell Type DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Lung" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix

se2.matrix <- as.matrix(experiment.se2.aggregate.filtered.normalized.dim.cluster@assays$RNA@counts)

es.max = sctype_score(scRNAseqData = se2.matrix, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(experiment.se2.aggregate.filtered.normalized.dim.cluster@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(experiment.se2.aggregate.filtered.normalized.dim.cluster@meta.data[experiment.se2.aggregate.filtered.normalized.dim.cluster@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(experiment.se2.aggregate.filtered.normalized.dim.cluster@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

# overlay cell types on UMAP

experiment.se2.aggregate.filtered.normalized.dim.cluster@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  experiment.se2.aggregate.filtered.normalized.dim.cluster@meta.data$customclassif[experiment.se2.aggregate.filtered.normalized.dim.cluster@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

png(file = "SE2 Cell Type Annotation.png", width = 600, height = 600) 
DimPlot(experiment.se2.aggregate.filtered.normalized.dim.cluster,
        reduction = "umap",
        label = TRUE,
        repel = TRUE,
        group.by = 'customclassif') + 
  plot_annotation(title = 'SE2 Cell Type Annotation')
dev.off()




#SE3

experiment.se3.aggregate.filtered.normalized.dim.cluster <- readRDS("~/data/R-scRNA-Seq-Processing/R Objects/experiment.se3.aggregate.filtered-normalized-dim-cluster.rds")

# get cell-type by cell matrix

se3.matrix <- as.matrix(experiment.se3.aggregate.filtered.normalized.dim.cluster@assays$RNA@counts)

es.max = sctype_score(scRNAseqData = se3.matrix, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(experiment.se3.aggregate.filtered.normalized.dim.cluster@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(experiment.se3.aggregate.filtered.normalized.dim.cluster@meta.data[experiment.se3.aggregate.filtered.normalized.dim.cluster@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(experiment.se3.aggregate.filtered.normalized.dim.cluster@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

# overlay cell types on UMAP

experiment.se3.aggregate.filtered.normalized.dim.cluster@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  experiment.se3.aggregate.filtered.normalized.dim.cluster@meta.data$customclassif[experiment.se3.aggregate.filtered.normalized.dim.cluster@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}


png(file = "SE3 Cell Type Annotation.png", width = 600, height = 600) 
DimPlot(experiment.se3.aggregate.filtered.normalized.dim.cluster,
        reduction = "umap",
        label = TRUE,
        repel = TRUE,
        group.by = 'customclassif') + 
  plot_annotation(title = 'SE3 Cell Type Annotation')
dev.off()


#all data
experiment.all.aggregate.filtered.normalized.dim.cluster <- readRDS("~/data/R-scRNA-Seq-Processing/R Objects/experiment.all.aggregate.filtered-normalized-dim-cluster.rds")




# get cell-type by cell matrix

all.matrix <- as.matrix(experiment.all.aggregate.filtered.normalized.dim.cluster@assays$RNA@counts)

es.max = sctype_score(scRNAseqData = all.matrix, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(experiment.all.aggregate.filtered.normalized.dim.cluster@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(experiment.all.aggregate.filtered.normalized.dim.cluster@meta.data[experiment.all.aggregate.filtered.normalized.dim.cluster@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(experiment.all.aggregate.filtered.normalized.dim.cluster@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

# overlay cell types on UMAP

experiment.all.aggregate.filtered.normalized.dim.cluster@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  experiment.all.aggregate.filtered.normalized.dim.cluster@meta.data$customclassif[experiment.all.aggregate.filtered.normalized.dim.cluster@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

png(file = "All Cell Type Annotation.png", width = 600, height = 600) 
DimPlot(experiment.all.aggregate.filtered.normalized.dim.cluster,
        reduction = "umap",
        label = TRUE,
        repel = TRUE,
        group.by = 'customclassif') + 
  plot_annotation(title = 'All Cell Type Annotation')
dev.off()




# Cell Type & Gene Interest Plots -----------------------------------------------
png(file = "Sftpc - SE3 - Cell Types.png", width = 600, height = 600)
VlnPlot(experiment.se3.aggregate.filtered.normalized.dim.cluster, 
        features = c("Sftpc"),
        split.by= "orig.ident",
        group.by = 'customclassif',
        pt.size=0)
dev.off()

png(file = "Lyz2 - SE3 - Cell Types.png", width = 600, height = 600)
VlnPlot(experiment.se3.aggregate.filtered.normalized.dim.cluster, 
        features = c("Lyz2"),
        split.by= "orig.ident",
        group.by = 'customclassif',
        pt.size=0)
dev.off()



png(file = "Sftpc - SE2 - Cell Types.png", width = 600, height = 600)
VlnPlot(experiment.se2.aggregate.filtered.normalized.dim.cluster, 
        features = c("Sftpc"),
        split.by= "orig.ident",
        group.by = 'customclassif',
        pt.size=0)
dev.off()

png(file = "Lyz2 - SE2 - Cell Types.png", width = 600, height = 600)
VlnPlot(experiment.se2.aggregate.filtered.normalized.dim.cluster, 
        features = c("Lyz2"),
        split.by= "orig.ident",
        group.by = 'customclassif',
        pt.size=0)
dev.off()




VlnPlot(experiment.se2.aggregate.filtered.normalized.dim.cluster, 
        features = c("Lyz2"),
        split.by= "orig.ident",
        group.by = 'customclassif',
        pt.size=0)


png(file = "Lyz2 - SE2 - Cell Types - Ridgeplot.png", width = 600, height = 600)
RidgePlot(experiment.se2.aggregate.filtered.normalized.dim.cluster,
          features = c("Lyz2"),
          group.by = 'customclassif')
dev.off()

png(file = "Lyz2 - SE3 - Cell Types - Ridgeplot.png", width = 600, height = 600)
RidgePlot(experiment.se3.aggregate.filtered.normalized.dim.cluster,
          features = c("Lyz2"),
          group.by = 'customclassif')
dev.off()


png(file = "Sftpc - SE2 - Cell Types - Ridgeplot.png", width = 600, height = 600)
RidgePlot(experiment.se2.aggregate.filtered.normalized.dim.cluster,
          features = c("Sftpc"),
          group.by = 'customclassif')
dev.off()

png(file = "Sftpc - SE3 - Cell Types - Ridgeplot.png", width = 600, height = 600)
RidgePlot(experiment.se3.aggregate.filtered.normalized.dim.cluster,
          features = c("Sftpc"),
          group.by = 'customclassif')
dev.off()




# Cell Nodes by Cluster -----------------------------------------------------


# load libraries
lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)


es.max = sctype_score(scRNAseqData = se3.matrix, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(experiment.se3.aggregate.filtered.normalized.dim.cluster@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(experiment.se3.aggregate.filtered.normalized.dim.cluster@meta.data[experiment.se3.aggregate.filtered.normalized.dim.cluster@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(experiment.se3.aggregate.filtered.normalized.dim.cluster@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])



# prepare edges
cL_resutls=cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

# prepare nodes
nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
for (i in 1:length(unique(cL_resutls$cluster))){
  dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

mygraph <- graph_from_data_frame(edges, vertices=nodes)

# Make the graph
gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")

scater::multiplot(DimPlot(experiment.se3.aggregate.filtered.normalized.dim.cluster, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss), gggr, cols = 2)





