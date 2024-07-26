library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(clusterProfiler)
library(org.Mm.eg.db)
library(igraph)
library(ggraph)
library(tidyverse)
library(fgsea)

########################
##PATHWAYS ANALYSIS
###Cell populations of interest:
#CD8_EM
#CD8_Naive
#CD8_Exhausted
#CD8_Effector
#CD8_CM1
#CD4_TH1_Helper
interested_populations <- c(
  "CD8_EM",
  "CD8_Naive",
  "CD8_Exhausted_2",
  "CD8_Effector",
  "CD8_CM1",
  "CD4_TH1_Helper"
)
#####Load the first dataframe:00119
df0019 <- read.table("Annotation_option_2_res_deg_0019.csv",header=T,sep=",")
df0034 <- read.table("Annotation_option_2_res_deg_0034.csv",header=T,sep=",")
df0037 <- read.table("Annotation_option_2_res_deg_0037.csv",header=T,sep=",")
df0119 <- read.table("Annotation_option_2_res_deg_0119.csv",header=T,sep=",")
df0307 <- read.table("Annotation_option_2_res_deg_0307.csv",header=T,sep=",")
dfpd1<- read.table("Annotation_option_2_res_deg_pd1.csv",header=T,sep=",")




#####FUNCTION FOR PATHWAYS ANALYSIS:::
####REACTOME::::
pathways_analysis_reactome <- function(dataframe, 
                                       cell_population_list, 
                                       upregulated_threshold, 
                                       downregulated_threshold, 
                                       treatment) {
  df <- data.frame(dataframe)
  combined_df <- data.frame()  # Initialize an empty data frame to store all results
  
  for (i in cell_population_list) {
    print(i)
    
    df_sub <- df %>% dplyr::filter(.id %in% as.character(i))
    
    ##################
    # Upregulated
    ##################
    df_sub_up <- df_sub[which(df_sub$avg_log2FC > as.numeric(upregulated_threshold) &
                                df_sub$p_val_adj < 0.01), ]
    df_sub_up_genes <- df_sub_up$Gene %>% as.vector() %>% as.character() %>% unique()
    ego_i_up <- enrichGO(gene = df_sub_up_genes, 
                         keyType = "SYMBOL", 
                         OrgDb = org.Mm.eg.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.1, 
                         readable = TRUE)
    ego_i_up_df <- data.frame(ego_i_up)
    if (nrow(ego_i_up_df) > 0) {
      ego_i_up_df <- ego_i_up_df[order(-ego_i_up_df$Count), ]
      ego_i_up_df$Regulation <- "UPREGULATED"
    }
    
    ##################
    # Downregulated
    ##################
    df_sub_down <- df_sub[which(df_sub$avg_log2FC < as.numeric(downregulated_threshold) &
                                  df_sub$p_val_adj < 0.01), ]
    df_sub_down_genes <- df_sub_down$Gene %>% as.vector() %>% as.character() %>% unique()
    ego_i_down <- enrichGO(gene = df_sub_down_genes, 
                           keyType = "SYMBOL", 
                           OrgDb = org.Mm.eg.db, 
                           ont = "BP", 
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.1, 
                           readable = TRUE)
    ego_i_down_df <- data.frame(ego_i_down)
    if (nrow(ego_i_down_df) > 0) {
      ego_i_down_df <- ego_i_down_df[order(-ego_i_down_df$Count), ]
      ego_i_down_df$Regulation <- "DOWNREGULATED"
    }
    
    ##################
    # Merge and Add Columns
    ##################
    merged <- rbind(ego_i_up_df, ego_i_down_df)
    if (nrow(merged) > 0) {
      merged$Cell_Population <- as.character(i)
      merged$Treatment <- as.character(treatment)
      
      ##################
      # Combine with Combined Data Frame
      ##################
      combined_df <- rbind(combined_df, merged)
    }
  }
  
  
  return(combined_df)
}
########################
####KEGG::::
########################
pathways_analysis_kegg <- function(dataframe, 
                                   cell_population_list, 
                                   upregulated_threshold, 
                                   downregulated_threshold, 
                                   treatment) {
  df <- data.frame(dataframe)
  combined_df <- data.frame()  # Initialize an empty data frame to store all results
  
  for (i in cell_population_list) {
    print(i)
    
    df_sub <- df %>% dplyr::filter(.id %in% as.character(i))
    
    ##################
    # Upregulated
    ##################
    df_sub_up <- df_sub[which(df_sub$avg_log2FC > as.numeric(upregulated_threshold) &
                                df_sub$p_val_adj < 0.01), ]
    df_sub_up_genes <- df_sub_up$Gene %>% as.vector() %>% as.character() %>% unique()
    print(df_sub_up_genes)
    df_sub_up_genes_entrez <- AnnotationDbi::select(org.Mm.eg.db, keys = df_sub_up_genes, 
                                                    columns = "ENTREZID", 
                                                    keytype = "SYMBOL")
    genes_up<-df_sub_up_genes_entrez$ENTREZID %>% as.vector() %>% na.omit() %>% unique()
    ego_i_up <- enrichKEGG(gene = genes_up,
                           organism = 'mmu',
                           pvalueCutoff = 0.5)
    ego_i_up_df <- data.frame(ego_i_up)
    if (nrow(ego_i_up_df) > 0) {
      ego_i_up_df <- ego_i_up_df[order(-ego_i_up_df$Count), ]
      ego_i_up_df$Regulation <- "UPREGULATED"
    }
    
    ##################
    # Downregulated
    ##################
    df_sub_down <- df_sub[which(df_sub$avg_log2FC < as.numeric(downregulated_threshold) &
                                  df_sub$p_val_adj < 0.01), ]
    df_sub_down_genes <- df_sub_down$Gene %>% as.vector() %>% as.character() %>% unique()
    print(df_sub_down_genes)
    df_sub_down_genes_entrez <- AnnotationDbi::select(org.Mm.eg.db, keys = df_sub_down_genes, 
                                                      columns = "ENTREZID", 
                                                      keytype = "SYMBOL")
    
    genes_down<-df_sub_down_genes_entrez$ENTREZID %>% as.vector() %>% na.omit() %>% unique()
    print(genes_down)
    ego_i_down <- enrichKEGG(gene = genes_down,
                             organism = 'mmu',
                             pvalueCutoff = 0.5)
    
    ego_i_down_df <- data.frame(ego_i_down)
    if (nrow(ego_i_down_df) > 0) {
      ego_i_down_df <- ego_i_down_df[order(-ego_i_down_df$Count), ]
      ego_i_down_df$Regulation <- "DOWNREGULATED"
    }
    
    ##################
    # Merge and Add Columns
    ##################
    merged <- rbind(ego_i_up_df, ego_i_down_df)
    if (nrow(merged) > 0) {
      merged$Cell_Population <- as.character(i)
      merged$Treatment <- as.character(treatment)
      
      ##################
      # Combine with Combined Data Frame
      ##################
      combined_df <- rbind(combined_df, merged)
      
    }
  }
  
  
  return(combined_df)
}


####write tables:
#####0119
df_0119_reactome<- pathways_analysis_reactome(df00119,interested_populations,1,-1,"0119") %>%
  mutate(Description = gsub(",", "|", Description))
name_r=paste0("0119","_","vs_PBS_REACTOME.csv")
write.table(df_0119_reactome,
            name_r,
            row.names=F,
            sep=",",
            quote=F)
df_0119_kegg<- pathways_analysis_kegg(df00119,interested_populations,0.5,-0.5,"0119") %>%
  mutate(Description = gsub(",", "|", Description))
name_k=paste0("0119","_","vs_PBS_KEGG.csv")
write.table(df_0119_kegg,
            name_k,
            row.names=F,
            sep=",",
            quote=F)

#####0034
df_0034_reactome<- pathways_analysis_reactome(df0034,interested_populations,1,-1,"0034") %>%
  mutate(Description = gsub(",", "|", Description))
name_r=paste0("0034","_","vs_PBS_REACTOME.csv")
write.table(df_0034_reactome,
            name_r,
            row.names=F,
            sep=",",
            quote=F)
df_0034_kegg<- pathways_analysis_kegg(df0034,interested_populations,0.5,-0.5,"0034") %>%
  mutate(Description = gsub(",", "|", Description))
name_k=paste0("0034","_","vs_PBS_KEGG.csv")
write.table(df_0034_kegg,
            name_k,
            row.names=F,
            sep=",",
            quote=F)




###0119
df_0119_reactome<- pathways_analysis_reactome(df0119,interested_populations,1,-1,"0119") %>%
  mutate(Description = gsub(",", "|", Description))
name_r=paste0("0119","_","vs_PBS_REACTOME.csv")
write.table(df_0119_reactome,
            name_r,
            row.names=F,
            sep=",",
            quote=F)
df_0119_kegg<- pathways_analysis_kegg(df0119,interested_populations,0.5,-0.5,"0119") %>%
  mutate(Description = gsub(",", "|", Description))
name_k=paste0("0119","_","vs_PBS_KEGG.csv")
write.table(df_0119_kegg,
            name_k,
            row.names=F,
            sep=",",
            quote=F)



##############################################################################################################################
####DATA VISUALIZATION
##############################################################################################################################
####PREPARING THE DATAFRAME::::
kegg_0019 <- read.table("0019_vs_PBS_KEGG.csv",header=T,sep=",")
kegg_0019$GeneRatio <- sapply(strsplit(kegg_0019$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
kegg_0019$Description <- sub(" -.*", "", kegg_0019$Description)

#####UPREGULATED
#####SELECT THE EM POPULATION:
kegg_0019_EM_up <- kegg_0019 %>% filter(Cell_Population %in% "CD8_EM") %>% filter(Regulation %in% "UPREGULATED") ##select upregulated and EM population
kegg_0019_EM_up <- kegg_0019_EM_up[order(-kegg_0019_EM_up$GeneRatio),] ##order by gene ratio
kegg_0019_EM_up_top10 <- kegg_0019_EM_up[1:10,] # select top ten
kegg_0019_EM_up_top10$Description <- sub(" -.*", "", kegg_0019_EM_up_top10$Description)

####DOT PLOT
jpeg("kegg_EM_up_top10.jpeg", width = 8, height = 6, units = "in", res = 300)
a <- ggplot(kegg_0019_EM_up_top10, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(
    title = "KEGG - Upregulated Genes in CD8-EM Population",
    x = "Gene Ratio",
    y = "Description",
    size = "Count",
    color = "p.adjust"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA)
  )
a
dev.off()

###BAR PLOT
jpeg("barplot_pathways_EM_UP.jpeg", width = 8, height = 6, units = "in", res = 300)
ggplot(kegg_0019_EM_up_top10, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(
    title = "KEGG - Upregulated Genes in CD8-EM Population",
    x = "Count",
    y = "Pathway Description",
    fill = "p.adjust"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA)
  )
dev.off()


######################
###NETWORK ANALYSIS
######################
####Generate a graph:
# Split geneID strings into individual genes
edges <- kegg_0019_EM_up_top10 %>%
  separate_rows(geneID, sep = "/") %>%
  select(Description, geneID)

# Create a list of unique pathways
unique_pathways <- unique(edges$Description)

# Create an adjacency matrix for pathways based on shared genes
adj_matrix <- matrix(0, nrow = length(unique_pathways), ncol = length(unique_pathways))
colnames(adj_matrix) <- rownames(adj_matrix) <- unique_pathways

for (gene in unique(edges$geneID)) {
  pathways_with_gene <- edges$Description[edges$geneID == gene]
  for (i in 1:length(pathways_with_gene)) {
    for (j in i:length(pathways_with_gene)) {
      if (i != j) {
        adj_matrix[pathways_with_gene[i], pathways_with_gene[j]] <- adj_matrix[pathways_with_gene[i], pathways_with_gene[j]] + 1
        adj_matrix[pathways_with_gene[j], pathways_with_gene[i]] <- adj_matrix[pathways_with_gene[j], pathways_with_gene[i]] + 1
      }
    }
  }
}

# Convert the adjacency matrix to a graph object
g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)

# Merge pathway attributes with graph nodes
V(g)$GeneRatio <- kegg_0019_EM_up_top10$GeneRatio[match(V(g)$name, kegg_0019_EM_up_top10$Description)]
V(g)$p.adjust <- kegg_0019_EM_up_top10$p.adjust[match(V(g)$name, kegg_0019_EM_up_top10$Description)]

# Plot the network graph
jpeg("enrichment_map.jpeg", width = 12, height = 8, units = "in", res = 300)
ggraph(g, layout = "fr") + 
  geom_edge_link(color = "grey") +
  geom_node_point(aes(size = GeneRatio, color = p.adjust)) +
  geom_node_text(aes(label = name), repel = TRUE) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(
    title = "Enrichment Map of Pathways",
    x = "",
    y = "",
    size = "Gene Ratio",
    color = "p.adjust"
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )
dev.off()


##########DOWNREGULATED:::::
#####SELECT THE EM POPULATION:
#####SELECT THE EM POPULATION:
kegg_0019_EM_down <- kegg_0019 %>% filter(Cell_Population %in% "CD8_EM") %>% filter(Regulation %in% "DOWNREGULATED") ##select downregulated and EM population
kegg_0019_EM_down <- kegg_0019_EM_down[order(-kegg_0019_EM_down$GeneRatio),] ##order by gene ratio
kegg_0019_EM_down_top10 <- kegg_0019_EM_down[1:10,] # select top ten
kegg_0019_EM_down_top10$Description <- sub(" -.*", "", kegg_0019_EM_down_top10$Description)

####DOT PLOT
jpeg("kegg_EM_down_top10.jpeg", width = 8, height = 6, units = "in", res = 300)
a <- ggplot(kegg_0019_EM_down_top10, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(
    title = "KEGG - downregulated Genes in CD8-EM Population",
    x = "Gene Ratio",
    y = "Description",
    size = "Count",
    color = "p.adjust"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA)
  )
a
dev.off()

###BAR PLOT
jpeg("barplot_pathways_EM_down.jpeg", width = 8, height = 6, units = "in", res = 300)
ggplot(kegg_0019_EM_down_top10, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(
    title = "KEGG - downregulated Genes in CD8-EM Population",
    x = "Count",
    y = "Pathway Description",
    fill = "p.adjust"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA)
  )
dev.off()


######################
###NETWORK ANALYSIS
######################
####Generate a graph:
# Split geneID strings into individual genes
edges <- kegg_0019_EM_down_top10 %>%
  separate_rows(geneID, sep = "/") %>%
  select(Description, geneID)

# Create a list of unique pathways
unique_pathways <- unique(edges$Description)

# Create an adjacency matrix for pathways based on shared genes
adj_matrix <- matrix(0, nrow = length(unique_pathways), ncol = length(unique_pathways))
colnames(adj_matrix) <- rownames(adj_matrix) <- unique_pathways

for (gene in unique(edges$geneID)) {
  pathways_with_gene <- edges$Description[edges$geneID == gene]
  for (i in 1:length(pathways_with_gene)) {
    for (j in i:length(pathways_with_gene)) {
      if (i != j) {
        adj_matrix[pathways_with_gene[i], pathways_with_gene[j]] <- adj_matrix[pathways_with_gene[i], pathways_with_gene[j]] + 1
        adj_matrix[pathways_with_gene[j], pathways_with_gene[i]] <- adj_matrix[pathways_with_gene[j], pathways_with_gene[i]] + 1
      }
    }
  }
}

# Convert the adjacency matrix to a graph object
g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)

# Merge pathway attributes with graph nodes
V(g)$GeneRatio <- kegg_0019_EM_down_top10$GeneRatio[match(V(g)$name, kegg_0019_EM_down_top10$Description)]
V(g)$p.adjust <- kegg_0019_EM_down_top10$p.adjust[match(V(g)$name, kegg_0019_EM_down_top10$Description)]

# Plot the network graph
jpeg("enrichment_map_down.jpeg", width = 12, height = 8, units = "in", res = 300)
ggraph(g, layout = "fr") + 
  geom_edge_link(color = "grey") +
  geom_node_point(aes(size = GeneRatio, color = p.adjust)) +
  geom_node_text(aes(label = name), repel = TRUE) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(
    title = "Enrichment Map of Pathways",
    x = "",
    y = "",
    size = "Gene Ratio",
    color = "p.adjust"
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )
dev.off()





##################################
####FIND COMMON PATHWAYS
##################################
find_common_descriptions <- function(df) {
  result <- list()
  
  # Unique Cell_Population values
  cell_populations <- unique(df$Cell_Population)
  
  for (cell_population in cell_populations) {
    # Subset data for the current Cell_Population
    subset_df <- df %>% filter(Cell_Population == cell_population)
    
    # Unique Regulation values
    regulations <- unique(subset_df$Regulation)
    
    # Initialize a list to store descriptions for each Regulation
    descriptions_list <- list()
    
    for (regulation in regulations) {
      # Subset data for the current Regulation
      reg_subset <- subset_df %>% filter(Regulation == regulation)
      
      # Store the descriptions
      descriptions_list[[regulation]] <- reg_subset$Description
    }
    
    # Find common descriptions across all Regulation groups
    common_descriptions <- Reduce(intersect, descriptions_list)
    
    # Store results
    result[[cell_population]] <- common_descriptions
  }
  
  return(result)
}


# Function to combine gene lists and sum numerical columns
combine_and_sum <- function(df) {
  result <- df %>%
    group_by(Description) %>%
    summarize(
      GeneRatio = mean(GeneRatio, na.rm = TRUE),
      BgRatio = unique(BgRatio),  # Assuming BgRatio should be unique
      pvalue = mean(pvalue, na.rm = TRUE),  # Taking minimum p-value
      p.adjust = mean(p.adjust, na.rm = TRUE),  # Taking minimum adjusted p-value
      qvalue = mean(qvalue, na.rm = TRUE),  # Taking minimum q-value
      geneID = paste(unique(unlist(strsplit(paste(geneID, collapse = "/"), "/"))), collapse = "/"),
      Count = sum(Count, na.rm = TRUE),
      Regulation = unique(Regulation),  # Assuming Regulation should be unique
      Cell_Population = unique(Cell_Population),  # Assuming Cell_Population should be unique
      Treatment = unique(Treatment),  # Assuming Treatment should be unique
      .groups = "drop"
    )
  return(result)
}


####RUN THE FUNCTION TO FIND OVERLAPPING GENE LISTS::
common_descriptions_per_cell_population <- find_common_descriptions(df=kegg_0019)

####AGAIN FOCUS ON THE EM POPULATION:::::
em_common_vector<- common_descriptions_per_cell_population[["CD8_EM"]] %>% as.vector()

#####SELECT THE EM POPULATION:
kegg_0019_EM_up <- kegg_0019 %>% filter(Cell_Population %in% "CD8_EM") ##select EM population
kegg_0019_EM_up <- kegg_0019_EM_up[order(-kegg_0019_EM_up$GeneRatio),] ##order by gene ratio

#####FILTER OVERLAPPING GENE SETS ONLY:
kegg_0019_EM_up_overlapping <- kegg_0019_EM_up %>% dplyr::filter(Description %in% em_common_vector)

###Combine overlapping terms:
kegg_0019_EM_up_overlapping_combined_df <- combine_and_sum(kegg_0019_EM_up_overlapping) %>% data.frame()
kegg_0019_EM_up_overlapping_combined_df_uniq <-kegg_0019_EM_up_overlapping_combined_df[duplicated(kegg_0019_EM_up_overlapping_combined_df$Description),]
kegg_0019_EM_up_overlapping_combined_df_uniq <- kegg_0019_EM_up_overlapping_combined_df_uniq[order(-kegg_0019_EM_up_overlapping_combined_df_uniq$Count),]

#####GENERATING BAR GRAPH AGAIN:

###BAR PLOT
jpeg("barplot_pathways_EM_overlapping.jpeg", width = 8, height = 6, units = "in", res = 300)
ggplot(kegg_0019_EM_up_overlapping_combined_df_uniq, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(
    title = "KEGG - Overlapping Pathways in CD8-EM Population",
    x = "Count",
    y = "Pathway Description",
    fill = "p.adjust"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA)
  )
dev.off()



#####PATHWAYS MAP:

######################
###NETWORK ANALYSIS
######################
####Generate a graph:
# Split geneID strings into individual genes
edges <- kegg_0019_EM_up_overlapping_combined_df_uniq %>%
  separate_rows(geneID, sep = "/") %>%
  select(Description, geneID)

# Create a list of unique pathways
unique_pathways <- unique(edges$Description)

# Create an adjacency matrix for pathways based on shared genes
adj_matrix <- matrix(0, nrow = length(unique_pathways), ncol = length(unique_pathways))
colnames(adj_matrix) <- rownames(adj_matrix) <- unique_pathways

for (gene in unique(edges$geneID)) {
  pathways_with_gene <- edges$Description[edges$geneID == gene]
  for (i in 1:length(pathways_with_gene)) {
    for (j in i:length(pathways_with_gene)) {
      if (i != j) {
        adj_matrix[pathways_with_gene[i], pathways_with_gene[j]] <- adj_matrix[pathways_with_gene[i], pathways_with_gene[j]] + 1
        adj_matrix[pathways_with_gene[j], pathways_with_gene[i]] <- adj_matrix[pathways_with_gene[j], pathways_with_gene[i]] + 1
      }
    }
  }
}

# Convert the adjacency matrix to a graph object
g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)

# Merge pathway attributes with graph nodes
V(g)$GeneRatio <- kegg_0019_EM_up_overlapping_combined_df_uniq$GeneRatio[match(V(g)$name, kegg_0019_EM_up_overlapping_combined_df_uniq$Description)]
V(g)$p.adjust <- kegg_0019_EM_up_overlapping_combined_df_uniq$p.adjust[match(V(g)$name, kegg_0019_EM_up_overlapping_combined_df_uniq$Description)]

# Plot the network graph
jpeg("enrichment_map_OVERLAPPING_lists.jpeg", width = 12, height = 8, units = "in", res = 300)
ggraph(g, layout = "fr") + 
  geom_edge_link(color = "grey") +
  geom_node_point(aes(size = GeneRatio, color = p.adjust)) +
  geom_node_text(aes(label = name), repel = TRUE) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(
       title = "Enrichment Map of Pathways",
       x = "",
       y = "",
       size = "Gene Ratio",
       color = "p.adjust"
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )
dev.off()



########################################################################################################################################################################################################################
#########TESTING
df00119_genes <- df00119[which(df00119$avg_log2FC > -1 &
                                 df00119$p_val_adj < 0.01), ]

genes <- df00119_genes$Gene 


df_sub_up_genes_entrez <- AnnotationDbi::select(org.Mm.eg.db, keys = genes, 
                                                columns = "ENTREZID", 
                                                keytype = "SYMBOL")
entrez<- df_sub_up_genes_entrez$ENTREZID %>% na.omit() %>% unique()
ego_down <- enrichKEGG(gene = entrez,
                       organism = 'mmu',
                       pvalueCutoff = 0.5)



######################
###ENRICHMENT PLOT
######################
kegg_0019_EM_up_top10$EnrichmentScore <- -log10(kegg_0019_EM_up_top10$p.adjust) * kegg_0019_EM_up_top10$Count
# Create a ranked list of genes
# Assuming 'Count' or another relevant column represents gene scores; adjust as needed
# Create a named vector of gene scores from the 'geneID' and 'Count' columns

# Create a data frame with gene scores
gene_scores <- kegg_0019_EM_up_top10 %>%
  rowwise() %>%
  mutate(genes = list(strsplit(as.character(geneID), "/")[[1]])) %>%
  unnest(genes) %>%
  ungroup() %>%
  mutate(score = Count) %>%
  group_by(genes) %>%
  summarize(score = sum(score))

# Convert data to a named vector
gene_scores_vector <- setNames(gene_scores$score, gene_scores$genes)

# Define gene sets for pathways
pathways <- list(
  "T cell receptor signaling pathway" = unlist(strsplit(kegg_0019_EM_up_top10$geneID[1], "/")),
  "Cell adhesion molecules" = unlist(strsplit(kegg_0019_EM_up_top10$geneID[2], "/")),
  "Cytokine-cytokine receptor interaction" = unlist(strsplit(kegg_0019_EM_up_top10$geneID[3], "/")),
  "Rheumatoid arthritis" = unlist(strsplit(kegg_0019_EM_up_top10$geneID[4], "/")),
  "Viral protein interaction with cytokine and cytokine receptor" = unlist(strsplit(kegg_0019_EM_up_top10$geneID[5], "/")),
  "Influenza A" = unlist(strsplit(kegg_0019_EM_up_top10$geneID[6], "/")),
  "Ribosome" = unlist(strsplit(kegg_0019_EM_up_top10$geneID[7], "/")),
  "Coronavirus disease" = unlist(strsplit(kegg_0019_EM_up_top10$geneID[8], "/")),
  "Malaria" = unlist(strsplit(kegg_0019_EM_up_top10$geneID[9], "/")),
  "Inflammatory bowel disease" = unlist(strsplit(kegg_0019_EM_up_top10$geneID[10], "/"))
)

# Perform GSEA analysis
fgseaRes <- fgsea(pathways = pathways, stats = gene_scores_vector,nperm=1000,scoreType = "pos")

# Generate Enrichment Score Plot for a specific pathway
jpeg("gsea_enrichment_plot.jpeg", width = 20, height = 8, units = "in", res = 300)
#plotEnrichment(fgseaRes, pathway = "T cell receptor signaling pathway")  # Adjust pathway name as needed
plotEnrichment(fgseaRes[["T cell receptor signaling pathway"]])
dev.off()

