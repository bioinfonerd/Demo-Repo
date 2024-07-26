

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




# Processing Normalization & Highest Variable Genes --------------------
location <- "~/data/R-scRNA-Seq-Processing/R_Objects/"

#Get List of R Objects to Process
rds_list <- list.files(path = location, pattern = "-fresh.rds", full.names = FALSE, recursive = FALSE)

for (i in rds_list){
  print(paste("Processing:",i))
  print("Loading Data")
  #load in data
  seurat.do <- readRDS(paste0(location,i))
  
  print("Normalizing")
  #normalize
  seurat.do <- SCTransform(seurat.do,
                           method = "glmGamPoi",
                           verbose = FALSE)
  
  print("Find Variable Features")
  #Find Variable Features
  seurat.do <- FindVariableFeatures(seurat.do, selection.method = "vst", nfeatures = 2000)
  
  #export Variable Features
  #VariableFeatures(seurat.do)
  
  # Identify the 10 most highly variable genes
  seurat.do.top10 <- head(VariableFeatures(seurat.do), 10)
  
  print("Plotting")
  # plot variable features with labels
  png(file = paste0("Variable Genes - ",
                    unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                    ".png"),
      width = 600, height = 600) 
  plot1 <- VariableFeaturePlot(seurat.do)
  print(LabelPoints(plot = plot1, points = seurat.do.top10, repel = TRUE))
  dev.off()
  
  print("Save R Object")
  #save
  saveRDS(seurat.do, file = paste0(unlist(strsplit(i,split="-"))[1],"-normalized.rds"))
}





