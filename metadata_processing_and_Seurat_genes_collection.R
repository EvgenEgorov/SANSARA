library(tidyverse)
library(Matrix)
library(Seurat)
saving_metadata <- function(pathway_to_metadata){
  dataset <- readRDS(pathway_to_metadata)
  dataset$barcode <- colnames(dataset)
  dataset$UMAP_1 <- dataset@reductions$umap@cell.embeddings[,1]
  dataset$UMAP_2 <- dataset@reductions$umap@cell.embeddings[,2]
  
  #Saving the metadata
  metadata <- dataset@meta.data
  write.csv(metadata, file='metadata.csv', quote=F, row.names=F)
  
  #Saving the count matrix
  counts_matrix <- GetAssayData(dataset, assay='RNA',slot='counts')
  writeMM(counts_matrix,file = 'counts.mtx')
  
  #Saving the PCA
  write.csv(dataset@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)
  
  #Saving gene names
  write.table(
    data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
    quote=F,row.names=F,col.names=F
  )
  
  #Saving variable features from Seurat
  write.csv(VariableFeatures(dataset),file='var_genes_seurat.csv',row.names = F)
}
