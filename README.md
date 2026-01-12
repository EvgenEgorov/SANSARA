# SANSARA — Splicing-Aware scrNa-Seq AppRoAch
SANSARA converts conventional scRNA-seq data into a splicing-aware gene expression matrix (saGEX), enabling analyses that account for gene-specific splicing dynamics at the single-cell level.

## Prerequisites

Before running SANSARA, ensure the following are installed:

- **Velocyto**: Follow installation instructions [here](https://velocyto.org/velocyto.py/install/index.html#install). Use Velocyto to generate loom files for your samples from Cell Ranger output as described in the [tutorial](https://velocyto.org/velocyto.py/tutorial/cli.html).  
- **Conda**: Create the SANSARA environment using:
  ```
  conda create --name <env> --file sansara_env.yml
  conda activate <env>
  ```
- **Additional Python packages**:
  ```
  pip install velovi torch
  ```
## Data Preparation
**Convert Seurat object to Anndata**  
   SANSARA requires counts, gene names, and metadata in Scanpy Anndata format. If needed you can extract these from a Seurat object using the provided R script:

```
Rscript pipeline/Seurat_data_extraction.R --input-file path/to/Seurat_rds_object 
```

## Running SANSARA
SANSARA uses `scvelo` and `velovi` to select informative genes and estimate RNA velocity.  

1. Edit the configuration in `pipeline/SANSARA_pipeline.py`:
   - Set file paths for counts, genes, metadata, and loom files
   - Ensure barcode naming is consistent
   - Adjust parameters if needed

2. Run the pipeline:
   ```
   python pipeline/SANSARA_pipeline.py
   ```
## Output
SANSARA generates splice-aware count matrix (saGEX), ready for downstream analyses such as normalization, clustering, dimensionality reduction, and differential expression testing.

## Citing

If you use SANSARA in your research, please cite the relevant publications:

- Daniil K Lukyanov, Evgeniy S Egorov, Valeriia V Kriukova, Denis Syrko, Victor V Kotliar, Kristin Ladell, David A Price, Andre Franke, Dmitry M Chudakov. "Splicing-aware scRNA-Seq resolution reveals execution-ready programs in effector Tregs". PLoS Comput Biol. 2025 
