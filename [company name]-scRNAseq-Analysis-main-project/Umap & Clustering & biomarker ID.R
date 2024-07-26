
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




# Processing UMAP, Clustering, and Biomarker Discovery --------------------
location <- "~/data/R-scRNA-Seq-Processing/R_Objects/"

#Get List of R Objects to Process
rds_list <- list.files(path = location, pattern = "-normalized.rds", full.names = FALSE, recursive = FALSE)


# UMAP and Cluster & Biomarker
for (i in rds_list){
  
  print(paste("Processing:",i))
  print("Loading Data")
  #load in data
  seurat.do <- readRDS(paste0(location,i))
  
  print("Dimension Reduction & Clustering")
  # Dimension Reduction & Clustering
  seurat.do <- RunPCA(seurat.do, verbose = FALSE)
  seurat.do <- RunUMAP(seurat.do, dims = 1:30, verbose = FALSE)
  seurat.do <- FindNeighbors(seurat.do, dims = 1:30, verbose = FALSE)
  seurat.do <- FindClusters(seurat.do, verbose = FALSE)

  #Plot UMAP + Clusters
  print("Plots")
  #All Samples
  png(file = paste0("Sample-Umap-Cluster-Overview-",
                    unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                    ".png"),
      width = 600, height = 600) 
  print(DimPlot(seurat.do,
          label = TRUE,
          raster = FALSE) + NoLegend())
  dev.off()
  
  
  #Plot by sample type UMAP + clustering
  png(file = paste0("Sample-Umap-Cluster-Overview-BySample",
                    unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                    ".png"),
      width = 600, height = 600) 
  print(DimPlot(seurat.do,
          label = TRUE,
          split.by= "orig.ident",
          raster = FALSE) + NoLegend())
  dev.off()
  
  
  # Find Markers (DONE) ------------------------------------------------------------

  # Only run on all dataset
  if(i == 'experiment.all_ids.filtered-normalized.rds'){
    print('Processing Marker Data for all Samples Comparison')
    for(i in unique(seurat.do@meta.data$orig.ident)){
      #all comparisons versus the controls
      output <- FindMarkers(seurat.do,ident.1 = "SE3_bleo",
                            ident.2 = i,
                            group.by = "orig.ident")
      
      #save table
      write.csv(output,file=paste0("Gene_Expression_SE3_bleo_vs_",i,".csv"))
      
      #all comparisons versus the controls
      output <- FindMarkers(seurat.do,ident.1 = "SE2_blemocyin",
                            ident.2 = i,
                            group.by = "orig.ident")
      
      #save table
      write.csv(output,file=paste0("Gene_Expression_SE2_bleo_vs_",i,".csv"))
      
      
      #all comparisons versus the controls
      output <- FindMarkers(seurat.do,ident.1 = "SE2_control",
                            ident.2 = i,
                            group.by = "orig.ident")
      
      #save table
      write.csv(output,file=paste0("Gene_Expression_SE2_ctrl_vs_",i,".csv"))
      
      #all comparisons versus the controls
      output <- FindMarkers(seurat.do,ident.1 = "SE3_control",
                            ident.2 = i,
                            group.by = "orig.ident")
      
      #save table
      write.csv(output,file=paste0("Gene_Expression_SE3_ctrl_vs_",i,".csv"))
    }
  }
  
  # Cluster Biomarkers -----------------------------------------------------
  
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  print("Cluster Biomarkers")
  seurat.do.markers <- FindAllMarkers(seurat.do,
                                      only.pos = TRUE,
                                      min.pct = 0.25,
                                      logfc.threshold = 0.25)
  
  
  #Export Markers As CSV
  write.csv(seurat.do.markers,
            file = paste0(unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                          ".csv"))


# Saving R Objects --------------------------------------------------------


  print("Saving R Objects")
  #save data
  saveRDS(seurat.do,
          file = paste0("experiment.",
                        unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                        ".aggregate.filtered.normalized.dim.cluster.rds"))
  
  #save data
  saveRDS(seurat.do.markers, 
          file = paste0("experiment.",
                        unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                        ".aggregate.filtered.normalized.dim.cluster.markers.rds"))
}

