from dataclasses import dataclass
from typing import Optional
import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re
import argparse
from pathlib import Path
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import torch
from velovi import preprocess_data, VELOVI

# ----------------------------
# Configuration container
# ----------------------------
@dataclass
class SansaraConfig:
    """
    Configuration object for running the SANSARA pipeline.

    Required:
    ---------
    loom_file : str
        Path to velocyto loom file.

    working_dir : str
        Directory containing counts.mtx, metadata.csv, gene_names.csv.

    Optional inputs:
    ----------------
    counts_file : str, default="counts.mtx"
        Name of count matrix file.

    metadata_file : str, default="metadata.csv"
        Name of metadata file.

    gene_names_file : str, default="gene_names.csv"
        Name of gene names file.

    barcode_prefix : Optional[str], default=None
        Prefix added to loom barcodes to match metadata.
        If None, prefix is inferred automatically.

    barcode_postfix : Optional[str], default=None
        Postfix added to loom barcodes to match metadata.
        If None, postfix is inferred automatically.

    output_file : str, default="splice_aware_matrix.csv"
        Output CSV file name.

    Analysis parameters:
    -------------------
    min_shared_counts : int, default=1
        Minimum shared counts between cells.

    n_top_genes : int, default=20000
        Number of genes selected for velocity analysis.

    n_pcs : int, default=30
        Number of PCA components.

    n_neighbors : int, default=30
        Number of neighbors for preprocessing.

    n_samples : int, default=25
        Number of samples used for velocity inference.

    latent_time_scaling : int, default=20
        Latent time scaling value for velocity inference.
    """

    loom_file: str
    working_dir: str

    counts_file: str = "counts.mtx"
    metadata_file: str = "metadata.csv"
    gene_names_file: str = "gene_names.csv"

    barcode_prefix: Optional[str] = None
    barcode_postfix: Optional[str] = None

    output_file: str = "splice_aware_matrix.csv"

    min_shared_counts: int = 1
    n_top_genes: int = 20000
    n_pcs: int = 30
    n_neighbors: int = 30
    n_samples: int = 25
    latent_time_scaling: int = 20

# ----------------------------
# Argparse
# ----------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "SANSARA converts conventional scRNA-seq data into a splicing-aware gene "
            "expression matrix (saGEX), enabling analyses that account for gene-specific "
            "splicing dynamics at the single-cell level."
        )
    )

    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument(
        "--loom-file", "-loom",
        required=True,
        help="Path to velocyto output file in loom format containing spliced and unspliced counts"
    )
    required_args.add_argument(
        "-wd", "--working-dir",
        required=True,
        help="Path to folder with counts, gene names, metadata (AnnData inputs)"
    )

    parser.add_argument("-counts", "--counts-file", default="counts.mtx")
    parser.add_argument("-metadata", "--metadata-file", default="metadata.csv")
    parser.add_argument("-genes", "--gene-names-file", default="gene_names.csv")

    parser.add_argument("-prefix", "--barcode-prefix", default=None)
    parser.add_argument("-postfix", "--barcode-postfix", default=None)

    parser.add_argument("-out", "--output-file", default="splice_aware_matrix.csv")

    parser.add_argument("-min-counts", "--min-shared-counts", default=1, type=int)
    parser.add_argument("-top", "--n-top-genes", default=20000, type=int)
    parser.add_argument("-pcs", "--n-pcs", default=30, type=int)
    parser.add_argument("-neighbors", "--n-neighbors", default=30, type=int)
    parser.add_argument("-samples", "--n-samples", default=25, type=int)
    parser.add_argument("-time", "--latent-time-scaling", default=20, type=int)

    return parser

def config_from_args(args: argparse.Namespace) -> SansaraConfig:
    return SansaraConfig(
        loom_file=args.loom_file,
        working_dir=args.working_dir,
        counts_file=args.counts_file,
        metadata_file=args.metadata_file,
        gene_names_file=args.gene_names_file,
        barcode_prefix=args.barcode_prefix,
        barcode_postfix=args.barcode_postfix,
        output_file=args.output_file,
        min_shared_counts=args.min_shared_counts,
        n_top_genes=args.n_top_genes,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        n_samples=args.n_samples,
        latent_time_scaling=args.latent_time_scaling,
    )

# ----------------------------
# Functions
# ----------------------------

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
    wd = Path(working_dir).expanduser().resolve()

    counts_path = wd / counts_file
    metadata_path = wd / metadata_file
    genes_path = wd / gene_names_file

    if not counts_path.exists():
        raise FileNotFoundError(counts_path)
    if not metadata_path.exists():
        raise FileNotFoundError(metadata_path)
    if not genes_path.exists():
        raise FileNotFoundError(genes_path)

    X = io.mmread(counts_path)
    adata = anndata.AnnData(X=X.transpose().tocsr())

    cell_meta = pd.read_csv(metadata_path)
    adata.obs = cell_meta
    adata.obs.index = adata.obs["barcode"]

    with open(genes_path, "r") as f:
        gene_names = f.read().splitlines()
    adata.var.index = gene_names

    return adata

def resolve_path(path_str: str, working_dir: str | None = None) -> Path:
    """
    Resolve path provided by used
    """
    p = Path(path_str).expanduser()
    if p.is_absolute():
        return p

    p_cwd = Path.cwd() / p
    if p_cwd.exists():
        return p_cwd.resolve()

    if working_dir is not None:
        p_wd = Path(working_dir).expanduser().resolve() / p
        return p_wd.resolve()

    return p_cwd.resolve()

def infer_prefix_postfix_from_adata(adata):
    """
    Infer prefix and postfix from adata.obs_names.
    Example:
        D01_rep1_AAACCTGAGGATTCGG-1
        -> prefix = "D01_rep1_"
        -> postfix = "-1"
    """
    s = str(adata.obs_names[0])

    # match "<prefix>_<barcode><optional -1>"
    m = re.match(r"^(.*_)([ACGTN]{10,})(-1)?$", s)

    prefix = m.group(1) if m else ""
    postfix = "-1" if s.endswith("-1") else ""

    return prefix, postfix

def process_loom_barcodes(loom_file, working_dir=None, prefix="", postfix=""):
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
    loom_path = resolve_path(loom_file, working_dir)
    if not loom_path.exists():
        raise FileNotFoundError(f"Loom file not found: {loom_path}")
    ldata = sc.read(str(loom_path), validate=False)

    def normalize(bc: str) -> str:
        bc = str(bc)
        if ":" in bc:
            bc = bc.split(":", 1)[1]
        bc = re.sub(r"x$", "", bc)  # remove trailing x
        return bc

    barcodes = [normalize(bc) for bc in ldata.obs_names]

    # apply prefix/postfix only if not already present
    if prefix or postfix:
        out = []
        for bc in barcodes:
            b = bc
            if prefix and not b.startswith(prefix):
                b = prefix + b
            if postfix and not b.endswith(postfix):
                b = b + postfix
            out.append(b)
        barcodes = out

    ldata.obs_names = barcodes
    ldata.var_names_make_unique()
    return ldata


def calculate_sagex_matrix(adata, velocity_layer="velocity"):
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
    genes_filtered_scvelo = adata.var.index[adata.var["velocity_genes"]].tolist()

    velocity = adata.to_df(layer=velocity_layer)[genes_filtered_scvelo]
    expression = adata.to_df()[genes_filtered_scvelo]

    velocity_expr_product = velocity * expression.values
    velocity_expr_product = velocity_expr_product.fillna(0)

    sagex_matrix = pd.DataFrame(index=velocity_expr_product.index)
    for gene in velocity_expr_product.columns:
        sagex_matrix[f"{gene}_unspliced"] = velocity_expr_product[gene].apply(
            lambda x: x if float(x) > 0 else 0
        )
        sagex_matrix[f"{gene}_spliced"] = velocity_expr_product[gene].apply(
            lambda x: (-1) * x if float(x) < 0 else 0
        )

    return sagex_matrix

# ----------------------------
# Analysis pipeline
# ----------------------------
def run_pipeline(cfg: SansaraConfig) -> str:
    """
    Run the full SANSARA splicing-aware analysis pipeline.

    Parameters
    ----------
    cfg : SansaraConfig
        Configuration object specifying input files, barcode mapping,
        preprocessing parameters, and output options.

        See `SansaraConfig` for the full list of available fields.

    Returns
    -------
    output_path : str
        Path to the generated splice-aware matrix CSV file.

    Example
    -------
    >>> cfg = SansaraConfig(
    ...     loom_file="demo.loom",
    ...     working_dir="./"
    ... )
    >>> out_path = run_pipeline(cfg)
    """
    print("Starting SANSARA analysis pipeline...")

    # Step 1
    print("\n[1/7] Loading count matrix and metadata...")
    adata = load_anndata(cfg.working_dir, cfg.counts_file, cfg.metadata_file, cfg.gene_names_file)
    print(f"  Loaded data: {adata.n_obs} cells × {adata.n_vars} genes")

    # Step 2
    print("\n[2/7] Loading loom file and processing barcode naming...")
    if cfg.barcode_prefix is not None or cfg.barcode_postfix is not None:
        prefix = cfg.barcode_prefix or ""
        postfix = cfg.barcode_postfix or ""
        print(f"Using user-provided barcode mapping: prefix='{prefix}', postfix='{postfix}'")
    else:
        prefix, postfix = infer_prefix_postfix_from_adata(adata)
        print(f"Auto-inferred barcode mapping: prefix='{prefix}', postfix='{postfix}'")
    ldata = process_loom_barcodes(
        cfg.loom_file,
        working_dir=cfg.working_dir,
        prefix=prefix,
        postfix=postfix,
    )
    print(f"  Loaded loom data: {ldata.n_obs} cells × {ldata.n_vars} genes")
    
    # Step 3
    print("\n[3/7] Merging AnnData and loom data...")
    # Ensure strings
    adata.obs_names = adata.obs_names.astype(str)
    ldata.obs_names = ldata.obs_names.astype(str)

    # Intersect cells first 
    common = adata.obs_names.intersection(ldata.obs_names)
    if len(common) == 0:
        raise ValueError(
            "No overlapping cell barcodes between counts/metadata and loom after processing. "
            "Printing adata.obs_names[:5] and ldata.obs_names[:5] to debug:"
        )
    adata = adata[common].copy()
    ldata = ldata[common].copy()

    # Merge
    adata = scv.utils.merge(adata, ldata)
    print(f"  Merged data: {adata.n_obs} cells × {adata.n_vars} genes")

    # Step 4
    print("\n[4/7] Preprocessing and gene selection...")
    scv.pp.filter_and_normalize(adata, min_shared_counts=cfg.min_shared_counts, n_top_genes=cfg.n_top_genes)
    sc.pp.pca(adata, n_comps=cfg.n_pcs)
    sc.pp.neighbors(adata, n_pcs=cfg.n_pcs, n_neighbors=cfg.n_neighbors)
    scv.pp.moments(adata)

    adata = preprocess_data(adata)

    # Step 5–6: VELOVI (make sure imports exist)
    print("\n[5/7] Training VeloVI model...")
    VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    vae = VELOVI(adata)
    vae.train()
    print("Model training complete")

    print("\n[6/7] Computing RNA velocities...")
    latent_time = vae.get_latent_time(n_samples=cfg.n_samples)
    velocities = vae.get_velocity(n_samples=cfg.n_samples, velo_statistic="mean")

    scaling = cfg.latent_time_scaling / latent_time.max(0)
    adata.layers["velocity"] = velocities / scaling

    # Step 7
    print("\n[7/7] Calculating splice-aware counts matrix...")
    sagex_matrix = calculate_sagex_matrix(adata)
    workdir = Path(cfg.working_dir).expanduser().resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    out_path = workdir / cfg.output_file
    sagex_matrix.to_csv(out_path, index=True)

    print("\n" + "=" * 70)
    print("Analysis complete!")
    print(f"saGEX matrix saved to: {out_path}")
    print(f"  Shape: {sagex_matrix.shape[0]} cells × {sagex_matrix.shape[1]} features")
    print(f"  ({len(sagex_matrix.columns)//2} genes with spliced/unspliced components)")
    print("=" * 70 + "\n")

    return out_path

# ----------------------------
# Execute pipeline 
# ----------------------------
def main(argv: Optional[list[str]] = None) -> str:
    parser = build_parser()
    args = parser.parse_args(argv)
    cfg = config_from_args(args)
    return run_pipeline(cfg)

if __name__ == "__main__":
    main()