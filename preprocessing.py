import os
import random

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy import sparse
from scipy.stats import median_abs_deviation

import plots


def is_outlier(adata: AnnData, metric: str, nmads: int = 5) -> pd.Series:
    """Identify outliers based on MAD (median absolute deviations).

    Parameters:
        adata (AnnData): Anndata object containing the data.
        metric (str): Name of the metric to check for outliers.
        nmads (int): Number of MADs to consider as threshold (default is 5).
    
    Returns:
        outlier (pd.Series): Boolean series indicating outliers.
        """
    M = adata.obs[metric]
    outlier = ((M < np.median(M) - nmads * median_abs_deviation(M)) |
               (np.median(M) + nmads * median_abs_deviation(M) < M))
    return outlier


def get_adata_from_multiple_samples(samples_path: str, n_samples: int, make_sparse: bool = False) -> AnnData:
    """Read h5ad matrices from a list of samples and concatenate them.
    
    Parameters:
        samples_path (str): Directory with h5ad matrices.
        n_samples (int): Number of samples to read.
        make_sparse (bool): Convert data to sparse matrix if True (default is False).
    
    Returns:
        adata (AnnData): Concatenated AnnData object from multiple samples.
    """
    samples_list = os.listdir(samples_path)
    mtx_list = random.sample(samples_list, k=n_samples)  # random sampling without replacement
    adatas_list = []

    for mtx in mtx_list:
        if os.path.isfile(os.path.join(samples_path, mtx)):
            adata = sc.read_h5ad(os.path.join(samples_path, mtx))
            # adata.obs_names_make_unique()
            if make_sparse:
                adata.X = sparse.csr_matrix(adata.X)
            adatas_list.append(adata)

    adata = ad.concat(adatas_list, join='outer')

    return adata

def lol():
    print('lol')

def qc_and_preprocess_to_pca(adata: AnnData, output_dir: str = None) -> AnnData:
    """
    Make QC and prerocessing of anndata matrix including PCA. 
    Args:
        adata (AnnData):
            Raw anndata matrix (obtained by combining matrices
            of individual donors and containing metadata).
        output_dir (str):
            If passed, save plots before and after QC to passed directory. 
    Returns:
        adata (AnnData):
            Processed anndata matrix.
        """
    # to fix the `TypeError` when writing to .h5ad 
    adata.obs['pool_number'] = [str(n) for n in adata.obs['pool_number']]
    adata.obs['age'] = [float(n) for n in adata.obs['age']]
    
    # 1. QC
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], inplace=True, percent_top=None, log1p=True)
    plots.plot_qc_distributions(adata, save_fig_to=output_dir, suffix='before')

    adata.obs["outlier"] = (is_outlier(adata, "log1p_total_counts") | is_outlier(adata, "log1p_n_genes_by_counts"))
    adata.obs['outlier'].value_counts()

    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3)
    adata.obs['mt_outlier'].value_counts()

    # filtering out:
    print(f"Total number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)]
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.filter_cells(adata, min_genes=200)
    print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

    plots.plot_qc_distributions(adata, save_fig_to=output_dir, suffix='after')

    # 2. PREPROCESSING
    adata.raw = adata.copy()
    adata.layers["counts"] = adata.X.copy()
    print('\nNormalizing, logarithming...')
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.layers["log_counts"] = adata.X.copy()

    

    print('\nSelect HVG...')
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=3000,
        flavor="seurat_v3",
        layer="counts"
    )

    adata = adata[:, adata.var.highly_variable]

    # sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"]) 
    sc.pp.scale(adata)
    sc.tl.pca(adata)

    return adata
