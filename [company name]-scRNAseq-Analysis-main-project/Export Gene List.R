gene_list_positive

cell_types <- c()
pos_gene_list <- c()

for(i in seq(1,length(gene_list_positive))){
  gene_set <- paste(gene_list_positive[[i]],collapse=",")
  pos_gene_list <- append(pos_gene_list, gene_set)
  cell_types <- append(cell_types, names(gene_list_positive[i]))
  
}

output <- data.frame(cell_types=cell_types,pos_gene_list=pos_gene_list)
write.csv(output,file="Cell_Type_Annotation_Definition.csv")
