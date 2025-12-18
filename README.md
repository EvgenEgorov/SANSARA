# SANSARA — Splicing-Aware scrNa-Seq AppRoAch
This page describes SANSARA pipeline to create splice-aware object from conventional scRNA-seq object. Before running SANSARA, velocyto software needs to be installed (https://velocyto.org/velocyto.py/install/index.html#install) and used on cellranger output for the samples of interest as described here https://velocyto.org/velocyto.py/tutorial/cli.html
SANSARA works in conda environment with the requirements listed in requirements.txt file. The environment can be created using the folowwing command:

```
conda create --name <env> --file requirements.txt 
```

After activating conda environment (see https://docs.conda.io/projects/conda/en/stable/user-guide/tasks/manage-environments.html), velovi and torch need to be installed via pip install in the environment. 

Conventional scRNAseq object (e.g. Seurat rds file) is converted to Scanpy Anndata from counts, gene names and metadata annotation in the pipeline. Optional R script to extract these data from Seurat object is provided in Seurat_data_extraction.R. Run 

```
Rscript pipeline/Seurat_data_extraction.R --input-file path/to/Seurat_rds_object 
```

SANSARA works with counts, genes and metadata tables and loom file from velocyto output. Pipeline uses scvelo and velovi to select informative set of genes and train the model to estimate velocity values. The output is saGEX count matrix, which can be used directly in the downstream single-cell analysis (normalisation, clustering, dimensional reduction, differential expression testing etc). To run SANSARA on your data, perform preparation steps described above, edit file paths, barcode naming (the naming should be aligned between expression data and loom file) and other parameters in pipeline/SANSARA_pipeline.py and run 

```
python pipeline/SANSARA_pipeline.py 
```
