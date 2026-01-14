## Install and load required packages
required_packages <- c("tidyverse", "Matrix", "Seurat", "argparser")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  install.packages(new_packages, repos = "https://cloud.r-project.org/")
}

suppressPackageStartupMessages(library(tidyverse))  
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(argparser))

## Parse arguments
p <- arg_parser("SANSARA - Export Seurat object to CSV/MTX format")
p <- add_argument(p, "--input_file", help="Path to Seurat object (in .rds format)", default=NA)
args <- parse_args(p)

if (is.na(args$input_file)) {
  stop("Error: --input_file is required")
}
if (!file.exists(args$input_file)) {
  stop(paste("Error: File not found:", args$input_file))
}

## Export data function
export_seurat_data <- function(seurat_path){
  message("Loading Seurat object...")
  dataset <- readRDS(seurat_path)
  dataset$barcode <- colnames(dataset)
  
  message("Extracting count matrix...")
  counts_matrix <- GetAssayData(dataset, assay='RNA', layer='counts')
  
  message("Saving metadata...")
  write.csv(dataset@meta.data, file='metadata.csv', quote=F, row.names=F)
  
  message("Saving count matrix...")
  writeMM(counts_matrix, file='counts.mtx')
  
  message("Saving gene names...")
  write.table(
    data.frame('gene'=rownames(counts_matrix)), file='gene_names.csv',
    quote=F, row.names=F, col.names=F
  )
   
  message(paste0("Successfully exported ", ncol(counts_matrix), " cells and ", 
                 nrow(counts_matrix), " genes"))
}

## Create output directory and export
original_dir <- getwd()
output_dir <- str_replace(str_extract(args$input_file, "[^/]*$"), "\\.[^.]*$", "")
dir.create(file.path(original_dir, output_dir), showWarnings = FALSE)
setwd(file.path(original_dir, output_dir))

export_seurat_data(args$input_file)

setwd(original_dir)
message(paste("Output saved to:", file.path(original_dir, output_dir)))
