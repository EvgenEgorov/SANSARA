import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import torch
from velovi import preprocess_data, VELOVI

# ============================================================================
# User parameters - change for your own data 
# ============================================================================

# Input paths
WORKING_DIR = 'path_to_sample/sample/' # Change this to your working directory
COUNTS_FILE = 'counts.mtx' 
METADATA_FILE = 'metadata.csv'
GENE_NAMES_FILE = 'gene_names.csv'
LOOM_FILE = 'path_to_loom/sample_alignments_G65EN.loom'

# Barcode naming configuration - the barcode names NEED to align between Anndata and loom files 
BARCODE_PREFIX = 'sample_'
BARCODE_POSTFIX = '-1'

# Output configuration
OUTPUT_FILE = 'splice_aware_matrix.csv'

# Analysis parameters
MIN_SHARED_COUNTS = 1
N_TOP_GENES = 20000
N_PCS = 30
N_NEIGHBORS = 30
N_SAMPLES = 25
LATENT_TIME_SCALING = 20

# ============================================================================
# Functions
# ============================================================================

def load_anndata(working_dir, counts_file, metadata_file, gene_names_file):
    """
    Load conventional count matrix, metadata, and gene names into AnnData object.
    
    Parameters:
    -----------
    working_dir : str
        Directory containing input files
    counts_file : str
        Path to counts matrix file (.mtx format)
    metadata_file : str
        Path to cell metadata file (.csv format)
    gene_names_file : str
        Path to gene names file (.csv format)
    
    Returns:
    --------
    adata : AnnData object
        
    """
    os.chdir(working_dir)
    
    # Load count matrix
    X = io.mmread(counts_file)
    adata = anndata.AnnData(X=X.transpose().tocsr())
    
    # Load and add metadata
    cell_meta = pd.read_csv(metadata_file)
    adata.obs = cell_meta
    adata.obs.index = adata.obs['barcode']
    
    # Load and add gene names
    with open(gene_names_file, 'r') as f:
        gene_names = f.read().splitlines()
    adata.var.index = gene_names
    
    return adata

def process_loom_barcodes(loom_file, prefix='', postfix=''):
    """
    Load loom file and adjust cell barcodes to match AnnData naming (change to your convention if needed).
    
    Parameters:
    -----------
    loom_file : str
        Path to loom file from velocyto output
    prefix : str
        Prefix to add to barcodes
    postfix : str
        Postfix to add to barcodes
    
    Returns:
    --------
    ldata : AnnData
        Loom data with adjusted barcodes
    """
    ldata = sc.read(loom_file, validate=False)
    
    # Extract and format barcodes
    barcodes = [bc.split(':')[1][:-1] for bc in ldata.obs.index]
    barcodes = [f"{prefix}{bc}{postfix}" for bc in barcodes]
    
    ldata.obs.index = barcodes
    ldata.var_names_make_unique()
    
    return ldata

def calculate_sagex_matrix(adata, velocity_layer='velocity'):
    """
    Calculate saGEX matrix.
    
    This function computes velocity * expression for each gene and splits
    the values into separate unspliced (positive) and spliced (negative) components.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object with velocity layer
    velocity_layer : str
        Name of the velocity layer in adata
    
    Returns:
    --------
    sagex_matrix : pd.DataFrame
        DataFrame with separate columns for spliced and unspliced components
    """
    # Get velocity genes
    genes_filtered_scvelo = adata.var.index[adata.var['velocity_genes']].tolist()
    
    # Extract velocity and expression data
    velocity = adata.to_df(layer=velocity_layer)[genes_filtered_scvelo]
    expression = adata.to_df()[genes_filtered_scvelo]
    
    # Calculate velocity * expression
    velocity_expr_product = velocity * expression.values
    velocity_expr_product = velocity_expr_product.fillna(0)
    
    # Split into unspliced (positive) and spliced (negative) components
    sagex_matrix = pd.DataFrame()
    for gene in velocity_expr_product.columns:
        # Unspliced: keep positive values
        sagex_matrix[f'{gene}_unspliced'] = velocity_expr_product[gene].apply(lambda x: x if float(x) > 0 else 0)
        # Spliced: convert negative values to positive
        sagex_matrix[f'{gene}_spliced'] = velocity_expr_product[gene].apply(lambda x: (-1)*x if float(x) < 0 else 0)
  
    return sagex_matrix

# ============================================================================
# Analysis pipeline
# ============================================================================

def main():
    """
    Main pipeline for splicing-aware count matrix calculation.
    """
    print("Starting SANSARA analysis pipeline...")
    
    # Step 1: Load AnnData object
    print("\n[1/7] Loading count matrix and metadata...")
    adata = load_anndata(WORKING_DIR, COUNTS_FILE, METADATA_FILE, GENE_NAMES_FILE)
    print(f"  Loaded data: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Step 2: Load and process loom data
    print("\n[2/7] Loading loom file and processing barcode naming...")
    ldata = process_loom_barcodes(LOOM_FILE, prefix=BARCODE_PREFIX, postfix=BARCODE_POSTFIX)
    print(f"  Loaded loom data: {ldata.n_obs} cells × {ldata.n_vars} genes")
    
    # Step 3: Merge AnnData object and loom data
    print("\n[3/7] Merging AnnData and loom data...")
    scv.utils.clean_obs_names(adata)
    scv.utils.clean_obs_names(ldata)
    adata = scv.utils.merge(adata, ldata)
    print(f"  Merged data: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Step 4: Preprocessing
    print("\n[4/7] Preprocessing and gene selection...")
    scv.pp.filter_and_normalize(adata, min_shared_counts=MIN_SHARED_COUNTS, n_top_genes=N_TOP_GENES)
    scv.pp.moments(adata, n_pcs=N_PCS, n_neighbors=N_NEIGHBORS)
    adata = preprocess_data(adata)
    
    # Step 5: Train VELOVI model
    print("\n[5/7] Training VeloVI model...")
    VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    vae = VELOVI(adata)
    vae.train()
    print("Model training complete")
    
    # Step 6: Calculate velocities
    print("\n[6/7] Computing RNA velocities...")
    latent_time = vae.get_latent_time(n_samples=N_SAMPLES)
    velocities = vae.get_velocity(n_samples=N_SAMPLES, velo_statistic="mean")
    
    # Scale velocities by latent time
    scaling = LATENT_TIME_SCALING / latent_time.max(0)
    adata.layers["velocity"] = velocities / scaling
    
    # Step 7: Calculate and save saGEX matrix
    print("\n[7/7] Calculating splice-aware conts matrix...")
    sagex_matrix = calculate_sagex_matrix(adata)
    sagex_matrix.to_csv(OUTPUT_FILE, index=True)
    
    print(f"\n{'='*70}")
    print(f"Analysis complete!")
    print(f"saGEX matrix saved to: {os.path.join(WORKING_DIR, OUTPUT_FILE)}")
    print(f"  Shape: {sagex_matrix.shape[0]} cells × {sagex_matrix.shape[1]} features")
    print(f"  ({len(sagex_matrix.columns)//2} genes with spliced/unspliced components)")
    print(f"{'='*70}\n")
    
#    return adata, sagex_matrix # optional

# ============================================================================
# Execute pipeline 
# ============================================================================

if __name__ == "__main__":
    main()
