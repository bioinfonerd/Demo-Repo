# Purpose: Creat Seurat Objects on [company name] scRNA-Seq Data


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


# Global Variables --------------------------------------------------------

#setup experiment conditions
experiment_name = "[company name] Example"
dataset_loc <- "/home/njohnson/data/aligned_data"
all_ids <- c("SE2_Drug_1","SE2_blemocyin","SE3_Drug_A",
             "SE3_bleo","_SE2_Drug_3A","SE2_Drug_2","SE2_control",
             "SE3_Drug_B","SE3_control")

se2_ids <- c("SE2_Drug_1","SE2_blemocyin","_SE2_Drug_3A","SE2_Drug_2","SE2_control")

se3_ids <- c("SE3_Drug_A","SE3_bleo","SE3_Drug_B","SE3_control")

se2_ctrl_vs_disease_ctrl <- c("SE2_blemocyin","SE2_control")

se3_ctrl_vs_disease_ctrl <- c("SE3_bleo","SE3_control")

se2_ctrl_vs_se3_ctrl <- c("SE2_control","SE3_control")

se2_disease_ctrl_vs_se3_disease_ctrl <- c("SE2_blemocyin","SE3_bleo")

b11_treatement_ctrl_disease_ctrl <- c("SE2_blemocyin","_SE2_Drug_3A","SE2_control")

se2_BCL_X_treatment_ctrl_disease_ctrl <- c("SE2_blemocyin","SE2_Drug_1","SE2_control")

se2_MCL_1_treatment_ctrl_disease_ctrl <- c("SE2_blemocyin","SE2_Drug_2","SE2_control")

se3_BCL_X_treatment_ctrl_disease_ctrl <- c("SE3_bleo","SE3_Drug_A","SE3_control")

se3_MCL_1_treatment_ctrl_disease_ctrl <- c("SE3_bleo","SE3_Drug_B","SE3_control")

all_BCL_X_treatment_ctrl_disease_ctrl <- c("SE2_blemocyin","SE2_Drug_1","SE2_control",
                                           "SE3_bleo","SE3_Drug_A","SE3_control")

all_MCL_1_treatment_ctrl_disease_ctrl <- c("SE2_blemocyin","SE2_Drug_2","SE2_control",
                                           "SE3_bleo","SE3_Drug_B","SE3_control")


# Experiment List 
Experiment_List <- list(all_ids,se2_ids,se3_ids,se2_ctrl_vs_disease_ctrl,
                     se3_ctrl_vs_disease_ctrl,se2_ctrl_vs_se3_ctrl,
                     se2_disease_ctrl_vs_se3_disease_ctrl,
                     b11_treatement_ctrl_disease_ctrl,
                     se2_BCL_X_treatment_ctrl_disease_ctrl,
                     se2_MCL_1_treatment_ctrl_disease_ctrl,
                     se3_BCL_X_treatment_ctrl_disease_ctrl,
                     se3_MCL_1_treatment_ctrl_disease_ctrl,
                     all_BCL_X_treatment_ctrl_disease_ctrl,
                     all_MCL_1_treatment_ctrl_disease_ctrl)

Experiment_Name_List <- c("all_ids","se2_ids","se3_ids","se2_ctrl_vs_disease_ctrl",
                     "se3_ctrl_vs_disease_ctrl","se2_ctrl_vs_se3_ctrl",
                     "se2_disease_ctrl_vs_se3_disease_ctrl",
                     "b11_treatement_ctrl_disease_ctrl",
                     "se2_BCL_X_treatment_ctrl_disease_ctrl",
                     "se2_MCL_1_treatment_ctrl_disease_ctrl",
                     "se3_BCL_X_treatment_ctrl_disease_ctrl",
                     "se3_MCL_1_treatment_ctrl_disease_ctrl",
                     "all_BCL_X_treatment_ctrl_disease_ctrl",
                     "all_MCL_1_treatment_ctrl_disease_ctrl")

names(Experiment_List) <- Experiment_Name_List

# Create Seurat Object ----------------------------------------------------

#load in each experiment data and all data to create separate Seurat objects 
# (uses alot of RAM, best to work with each object and save-current instance is large enough so won't separate)

#load in data
for (i in seq(1,length(Experiment_List))){
  print(names(Experiment_List[i]))
  print(i)
  
  #loading in data
  print("Loading In Data")
  bc.matrix.data <- lapply(Experiment_List[[i]], function(i){
    d10x <- Read10X_h5(file.path(dataset_loc,paste0(i,"/outs"),"raw_feature_bc_matrix.h5"))
    colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
    d10x
  })
  names(bc.matrix.data) <- names(Experiment_List[i])
 
  #combine dataset
  print("Combine Data")
  all.data <- do.call("cbind", bc.matrix.data)
  
  
  #create Seurat Object
  print("Creating Seurat Object")
  #filter criteria: 
  #remove genes that do not occur in a minimum of 0 cells
  #remove cells that don’t have a minimum of 200 features
  seurat.object <- CreateSeuratObject(
    all.data,
    min.cells = 1,
    min.features = 300,
    names.field = 2,
    names.delim = "\\-")
  
  #MIT genes
  print("Assigning Mitochondrial Status")
  seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, pattern = "^MT-")
  
  #save object
  print("Saving Seurat Object")
  saveRDS(seurat.object, file = paste0("experiment.",names(Experiment_List[i]),".filtered-fresh.rds"))
  
}
