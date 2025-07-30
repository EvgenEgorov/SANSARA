suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(argparser))

p <- arg_parser("SANSARA")
p <- add_argument(p,"--input_file",help="Seurat object (in .rds format)",default=NA)
args <- parse_args(p)

saving_metadata <- function(pathway_to_metadata){
  dataset <- readRDS(pathway_to_metadata)
  dataset$barcode <- colnames(dataset)

#Saving the metadata
metadata <- dataset@meta.data
write.csv(metadata, file='metadata.csv', quote=F, row.names=F)

#Saving the count matrix
counts_matrix <- GetAssayData(dataset, assay='RNA',layer='counts')
writeMM(counts_matrix,file = 'counts.mtx')

#Saving the PCA
#  write.csv(dataset@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)

#Saving gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)

#Saving variable features from Seurat
#  write.csv(VariableFeatures(dataset),file='var_genes_seurat.csv',row.names = F)
}
directory <- getwd()
output_dir <- str_replace(str_extract(args$input_file, "[^/]*$"), "\\.[^.]*$", "")
dir.create(file.path(directory,output_dir), showWarnings = FALSE)
setwd(file.path(directory,output_dir))
saving_metadata(args$input_file)
