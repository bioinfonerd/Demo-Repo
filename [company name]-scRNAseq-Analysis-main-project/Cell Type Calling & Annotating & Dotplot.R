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






# Processing Cell Type Calling & Annotation --------------------

location <- "~/data/R-scRNA-Seq-Processing/R_Objects/"

#completed
#"experiment.all.aggregate.filtered.normalized.dim.cluster.rds"
#"experiment.se2.aggregate.filtered.normalized.dim.cluster.rds"

#Get List of R Objects to Process
rds_list <- list.files(path = location, pattern = ".dim.cluster.rds", full.names = FALSE, recursive = FALSE)

# Cell Type DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Lung" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)


#provided gene set
Kobayashi_Senescence_Signature <- c("Calm1","Calm2","Capn2","Ccnd1","Ccnd2","Cdk4","Cdkn1a","Cdkn2b","Eif4ebp1",
  "Ets1","Gadd45b","H2-D1","H2-K1","h2-Q4","H2-Q6","H2-Q7","H2-T22","H2-T23",
  "Kras","Map2k2","Mapk13","Mapk3","Myc","Ppid","Ppp1ca","Ppp3ca","Rbbp4",
  "Rheb","Rras","Sqstm1","Slc25a4","Tgfb2","Trp53","Vdac1","Vdac2","Vdac3",
  "Zfp36l1")

Simple_Senescence_Signature <- c("p16","p15","p19","p21","p27","PAI-1")

#make list of a list for cell type definition
gs_provided_list_positive = list(Kobayashi_Senescence_Signature,Simple_Senescence_Signature)
names(gs_provided_list_positive) <- c("Kobayashi_Senescence_Signature","Simple_Senescence_Signature")
gs_provided_list_negative <- list(c(""),c(""))
names(gs_provided_list_negative) <- c("Kobayashi_Senescence_Signature","Simple_Senescence_Signature")

gs_provided_list <- list(gs_provided_list_positive,gs_provided_list_negative)
names(gs_provided_list) <- c("gs_positive","gs_negative")

# Gene Sets of Interest (Dot Plot)
Profibrotic_Genes <- c("Col1a1", "Acta2", "Timp1", "Serpine1", "Fn1")
Senescence_related_genes <- c("Cdkn1a", "Cdkn1b", "Cdkn2b",
                              "Glb1","Serpine1")
SASP_genes <- c("Il6", "Ccl2", "Fgf2","Tgfb1",
                "Egf","Cxcl1","Col1a1")
Ephrin_Ligands <- c("Efna1", "Efna2", "Efna3","Efna4", "Efna5","Efnb1",
                    "Efnb2","Efnb3")
Ephrin_Receptor <- c("Epha1", "Epha2", "Epha3","Epha4", "Epha5","Epha6",
                     "Epha7","Epha8","Epha10","Ephb1","Ephb2","Ephb3",
                     "Ephb4","Ephb6")
BCL2_protein_family <- c("Bax", "Bak1", "Bok",
                         "Bcl2l11", "Bid","Bcl2",
                         "Bcl2l2","Bcl2l1",
                         "Mcl1","Bcl2l10","Bbc3",
                         "Bmf","Bad","Pmaip1","Hrk",
                         "Bnip3","Bik")
Condition_474 <- c("Nr1d1", "Abtb1",
                   "Dbp", "Prdm8",
                   "1110038F14Rik","Nfil3","H2-Q6",
                   "Sppl3","Bhlhe41","Nr4a3",
                   "Ccdc85b","Pip4p1","Sema6c",
                   "Gm8229","Gata3","Mcpt4","Tmem252")
Condition_124 <- c("Rbm25","Msl3l2",
                   "Rnf187", "Stub1","Nr1d1",
                   "Nub1","Cops8",
                   "Tmem41a","Mtch1",
                   "Thap11","Kif20b","Chpf",
                   "Igkv3-2","Ppp2cb","Emc10","Zfp574")
Condition_394 <- c("Nr1d1", "Tmem252", "Camkk1",
                   "Plekhf1","Id2",
                   "Nr4a3","Msl3l2","Bhlhe41",
                   "Adm","Nfil3","Cebpd",
                   "Kctd11","Socs3","Tmem8","Slc38a2",
                   "1200007C13Rik","Dbp","Myc","4930516B21Rik")
Condition_470 <- c("Arntl", "Dbp", "Irf2bp2",
                   "Hspa5","Bhlhe41",
                   "Npas2","Tmem252","Pip4p1",
                   "Hspa8","Adm","Nfil3",
                   "Stub1","Ywhaz","Rrbp1","Sec31a",
                   "Tomm20","Rasl11a","Id2","Rab17")

# Gene List 
gene_list <- c(Profibrotic_Genes,Senescence_related_genes,SASP_genes,
               Ephrin_Ligands,Ephrin_Receptor,BCL2_protein_family,
               Condition_474,Condition_124,Condition_394,Condition_470)

# Genes of Interest
goi <- c("Sftpc","Lyz2","Sftpa1")

# For Loop Processing -------------------------------------------

for (i in rds_list){

  print(paste("Processing:",i))
  print("Loading Data")
  #load in data
  seurat.do <- readRDS(paste0(location,i))

# Cell Type Annotation (ScType) (DONE) -----------------------------------------
  print('Cell Type Annotation')
  #https://www.nature.com/articles/s41467-022-28803-w, https://github.com/IanevskiAleksandr/sc-type/
  # get cell-type by cell matrix
  seurat.matrix <- as.matrix(seurat.do@assays$RNA@counts)
  
  es.max = sctype_score(scRNAseqData = seurat.matrix, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
  
  # merge by cluster
  cL_resutls = do.call("rbind", lapply(unique(seurat.do@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seurat.do@meta.data[seurat.do@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat.do@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  #print(sctype_scores[,1:3])
  
  # overlay cell types on UMAP
  seurat.do@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    seurat.do@meta.data$customclassif[seurat.do@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }

# Cell Type Annotation Plot (DONE) -----------------------------------------------
  print('Cell Type Annotation')
  png(file = paste0("Cell Type_Annotation_Umap-",
                    unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                    ".png"),
      width = 600, height = 600) 
  print(DimPlot(seurat.do,
          reduction = "umap",
          label = TRUE,
          repel = TRUE,
          group.by = 'customclassif') + 
    ggtitle(unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2]) +
    theme(plot.title = element_text(hjust = 0.5)))
  dev.off()

# Cell Type by Cluster & Sample (DONE)------------------------------------------
  print('Cluster & Sample Plots')
  #group by sample id, clusters, cell type
  clusters.cell_type <- seurat.do@meta.data %>% 
    group_by(orig.ident,seurat_clusters, customclassif) %>% 
    summarise(n = n())
  
  #plot
  png(file = paste0("Cell Type_Cluster_",
                    unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                    ".png"),
      width = 600, height = 600)
  
  print(ggplot(clusters.cell_type, aes(x=seurat_clusters,
                                 y=n, 
                                 fill=orig.ident)) + 
    geom_bar(position="dodge", stat="identity") +
    facet_wrap(~customclassif) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
    scale_fill_viridis(discrete = T, option = "E"))
  
  dev.off()

# Dotplot (Gene List) (DONE) ---------------------------------------------------
  print('Print Dotplot')
  gene_sets <- list(Profibrotic_Genes,Senescence_related_genes,SASP_genes,
                    Ephrin_Ligands,Ephrin_Receptor,BCL2_protein_family,
                    Condition_474,Condition_124,Condition_394,Condition_470)
  
  title <- list("Profibrotic Genes","Senescence Related Genes","SASP Genes",
                "Ephrin Ligands","Ephrin Receptor","BCL2 protein family",
                "Condition 474","Condition 124","Condition 394",
                "Condition 470")
  
  for (n in seq(from=1, to=length(gene_sets))){
    png(file = paste0(title[[n]],
                      "-Dotplot-",
                      unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                      ".png"),
        width = 600, height = 600)
    
    print(DotPlot(seurat.do,
                  features = gene_sets[[n]],
                  cols = c("blue", "red"),
                  group.by = "orig.ident") + 
            RotatedAxis() +
            ggtitle(title[[n]]) +
            theme(plot.title = element_text(hjust = 0.5)))
    dev.off()
    
  }
  

  
# Gene of Interest Plots (GOI) (DONE)-------------------------------------------
  print("Gene of Interest")
  for (n in goi){
    
    #fails if not found, check then proceed if not
    if (n %in% rownames(seurat.do@assays$RNA@counts)) {
      print('Gene is Expressed')
    } else {
      print('Gene is NOT Expressed')
      next
    }
    
    # Violin Plot of Gene Across Clusters
    png(file = paste0("VlnPlot-Expression across Clusters-GOI-",
                      n,
                      "-",
                      unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                      ".png"), width = 600, height = 600) 
    print(VlnPlot(seurat.do,
            features = n,
            pt.size=0,
            ncol = 1)) 
    dev.off()
    
    # Violin Plot of Gene Across Clusters
    png(file = paste0("VlnPlot-Expression across Clusters & Samples-GOI-",
                      n,
                      "-",
                      unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                      ".png"), width = 600, height = 600) 
    print(VlnPlot(seurat.do,
            features = n,
            split.by= "orig.ident",
            pt.size=0,
            ncol = 1))
    dev.off()

    # Violin Plot Across Cell Types & Samples
    png(file = paste0("VlnPlot-Expression across Cell Types & Samples-GOI-",
                      n,
                      "-",
                      unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                      ".png"), width = 600, height = 600) 
    print(VlnPlot(seurat.do, 
            features = n,
            split.by= "orig.ident",
            group.by = 'customclassif',
            pt.size=0))
    dev.off()
    
    # Gene Heat Map across UMAP
    png(file = paste0("HeatMap-UMAP-GOI-",
                      n,
                      "-",
                      unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                      ".png"), width = 600, height = 600) 
    print(FeaturePlot(seurat.do,
                features = n,
                pt.size = 0.2,
                ncol = 1))
    dev.off()

    # Ridge Plot Across Cell Types & All Samples
    png(file = paste0("Ridgeplot-Cell Types-All Samples-GOI-",
                      n,
                      "-",
                      unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                      ".png"), width = 600, height = 600) 
    print(RidgePlot(seurat.do,
              features = n,
              group.by = 'customclassif'))
    dev.off()
  }
  
  
  

# Gene Lists Plots (DONE)-------------------------------------------------------
  print("Gene List Plots")

  
  for (n in gene_list){
    
      #fails if not found, check then proceed if not
  if (n %in% rownames(seurat.do@assays$RNA@counts)) {
    print('Gene is Expressed')
  } else {
    print('Gene is NOT Expressed')
    next
  }
    
    
    # Violin Plot of Gene Across Clusters
    png(file = paste0("VlnPlot-Expression across Clusters-Gene List-",
                      n,
                      "-",
                      unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                      ".png"), width = 600, height = 600) 
    print(VlnPlot(seurat.do,
                  features = n,
                  pt.size=0,
                  ncol = 1)) 
    dev.off()
    
    # Violin Plot of Gene Across Clusters
    png(file = paste0("VlnPlot-Expression across Clusters & Samples-Gene List-",
                      n,
                      "-",
                      unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                      ".png"), width = 600, height = 600) 
    print(VlnPlot(seurat.do,
                  features = n,
                  split.by= "orig.ident",
                  pt.size=0,
                  ncol = 1))
    dev.off()
    
    # Violin Plot Across Cell Types & Samples
    png(file = paste0("VlnPlot-Expression across Cell Types & Samples-Gene List-",
                      n,
                      "-",
                      unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                      ".png"), width = 600, height = 600) 
    print(VlnPlot(seurat.do, 
                  features = n,
                  split.by= "orig.ident",
                  group.by = 'customclassif',
                  pt.size=0))
    dev.off()
    
    # Gene Heat Map across UMAP
    png(file = paste0("HeatMap-UMAP-Gene List-",
                      n,
                      "-",
                      unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                      ".png"), width = 600, height = 600) 
    print(FeaturePlot(seurat.do,
                      features = n,
                      pt.size = 0.2,
                      ncol = 1))
    dev.off()
    
    # Ridge Plot Across Cell Types & All Samples
    png(file = paste0("Ridgeplot-Cell Types-All Samples-Gene List-",
                      n,
                      "-",
                      unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                      ".png"), width = 600, height = 600) 
    print(RidgePlot(seurat.do,
                    features = n,
                    group.by = 'customclassif'))
    dev.off()
  }


# Save R Object (DONE)----------------------------------------------------------
  print("Save R Object")
  saveRDS(seurat.do,
        file = paste0("experiment.",
                      unlist(strsplit(unlist(strsplit(i,split="-"))[1],split="[.]"))[2],
                      ".aggregate.filtered.normalized.dim.cluster.celltype.rds"))
}
  
# Cell Type Calling By Provided Cell Type Definitions List --------------------------------------
  # 
  # es.max = sctype_score(scRNAseqData = seurat.matrix, scaled = TRUE, 
  #                       gs = gs_provided_list$gs_positive, gs2 = gs_provided_list$gs_negative) 
  # 
  # # merge by cluster
  # cL_resutls = do.call("rbind", lapply(unique(seurat.do@meta.data$seurat_clusters), function(cl){
  #   es.max.cl = sort(rowSums(es.max[ ,rownames(seurat.do@meta.data[seurat.do@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  #   head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat.do@meta.data$seurat_clusters==cl)), 10)
  # }))
  # sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  # 
  # # set low-confident (low ScType score) clusters to "unknown"
  # sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  # print(sctype_scores[,1:3])
  # 
  # # overlay cell types on UMAP
  # seurat.do@meta.data$customclassif = ""
  # for(j in unique(sctype_scores$cluster)){
  #   cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  #   seurat.do@meta.data$customclassif[seurat.do@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  # }
  



  



