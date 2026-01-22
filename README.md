# SANSARA — Splicing-Aware scrNa-Seq AppRoAch
SANSARA converts conventional scRNA-seq data into a splicing-aware gene expression matrix (saGEX), enabling analyses that account for gene-specific splicing dynamics at the single-cell level.

## Installation
  
Clone the repositoty to access SANSARA scripts and demo files:
  ```
  git clone https://github.com/EvgenEgorov/SANSARA.git
  cd SANSARA
  ``` 
Create SANSARA environment with Conda:
  ```
  conda env create --name sansara_env --file sansara_env.yml
  conda activate sansara_env
  ```
Check installation by running help command:
  ```
  python pipeline/SANSARA_pipeline.py --help
  ```
## Data Preparation

Before running SANSARA, ensure the following files exist:

**Velocyto loom file**: Follow installation instructions [here](https://velocyto.org/velocyto.py/install/index.html#install). Use Velocyto to generate loom files for your samples from Cell Ranger output as described in the [tutorial](https://velocyto.org/velocyto.py/tutorial/cli.html).

**Convert Seurat object to Anndata**  
SANSARA requires counts, gene names, and metadata in Scanpy Anndata format. If needed you can extract these from a Seurat object using the provided R script:

```
Rscript pipeline/Seurat_data_extraction.R --input_file path/to/Seurat_rds_object 
```
R script will create a folder named after Seurat .rds object. It should contain matrix.mtx, gene\_names.csv and metadata.csv files required for running SANSARA. 

## Running SANSARA
SANSARA uses `scvelo` and `velovi` to select informative genes and estimate RNA velocity.  
The pipeline can be executed in two ways:  
  - Option 1 — Command line interface  
    Best for batch runs, pipelines, and reproducible execution.  
  - Option 2 — Jupyter Notebook   
    Best for interactive exploration, debugging, and parameter tuning.      
  
### Option 1 — Command line interface  
**Basic usage**
```
python pipeline/SANSARA_pipeline.py --loom-file path/to/object.loom -wd path/to/anndata_folder
```  
**Required arguments**
    
`--loom-file <file.loom>`   
    Path to the output of velocyto in .loom format  
`-wd <path_adata> `  
    Path to the folder with counts, gene names, and metadata in Scanpy Anndata format  

**Optional arguments**
- Input files
```
--counts-file        Name of count matrix file, default = 'counts.mtx'
--metadata-file      Name of metadata file, default = 'metadata.csv'
--gene-names-file    Name of gene names file, default = 'counts.mtx'
```
- Barcode naming: barcode names should align between Anndata and loom files
```
--barcode-prefix     Prefix added to loom barcodes to match metadata. If None, prefix is inferred automatically.
--barcode-postfix    Postfix added to loom barcodes to match metadata. If None, postfix is inferred automatically.
```   
- Output
```
--output-file        Output CSV file name, default = "splice_aware_matrix.csv"
```
- Analysis parameters
```
--min-shared-counts     Minimum shared counts between cells in the sample, default = 1
--n-top-genes           Number of genes selected for the velocity analysis, default = 20000 
--n-pcs                 Number of principal components used for preprocessing, default = 30
--n-neighbors           Number of neighbors used for in the preprocessing, default = 30
--n-samples             Number of samples used for velocity inference, default = 25
--latent-time-scaling   Latent time scaling value used in the velocity inference, default = 20
```
### Option 2 — Run from Jupyter Notebook  
SANSARA can also be executed interactively from a notebook.  
Example SANSARA_demo.ipynb can be found in `demo` folder.  
**Example**
```
from pipeline import SANSARA_pipeline as sansara

cfg = sansara.SansaraConfig(
    loom_file="path/to/object.loom",
    working_dir="path/to/anndata_folder",
)

output_path = sansara.run_pipeline(cfg)
```

## Output
SANSARA generates a splice-aware gene expression matrix (saGEX) saved as a CSV file.  
The output is ready for downstream analyses such as normalization, clustering, dimensionality reduction, and differential expression testing.

## Citing

If you use SANSARA in your research, please cite the relevant publications:

- Daniil K Lukyanov, Evgeniy S Egorov, Valeriia V Kriukova, Denis Syrko, Victor V Kotliar, Kristin Ladell, David A Price, Andre Franke, Dmitry M Chudakov. "Splicing-aware scRNA-Seq resolution reveals execution-ready programs in effector Tregs". PLoS Comput Biol.2025  
[doi:10.1371/journal.pcbi.1013682](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013682)
