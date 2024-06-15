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


def aggregate_and_filter(
    adata,
    cell_identity=None,
    donor_key="bulk_sample",
    condition_key="dea_labels",
    cell_identity_key=None,
    obs_to_keep=[],  # which additional metadata to keep, e.g. gender, age, etc.
    replicates_per_patient=1,
    min_cells_per_donor=30
    ):
    """
    Creates an AnnData object with one pseudo-replicate for each donor for a specified 
    subpopulation from the original single-cell AnnData object. 
    Also it filters out donors that have fewer than `min_cells_per_donor' cells 
    for the specified population.

    By changing the `replicates_per_patient` parameter, 
    several (n) pseudo-replicates can be created for each sample; 
    cells are then split into n subsets of roughly equal sizes.
    """
    # subset adata to the given cell identity
    if cell_identity:
        adata_cell_pop = adata[adata.obs[cell_identity_key] == cell_identity].copy()
    adata_cell_pop = adata
    
    # check which donors to keep according to the number of cells specified with NUM_OF_CELL_PER_DONOR
    size_by_donor = adata_cell_pop.obs.groupby([donor_key]).size()
    donors_to_drop = [
        donor
        for donor in size_by_donor.index
        if size_by_donor[donor] <= min_cells_per_donor
    ]
    if len(donors_to_drop) > 0:
        print("Dropping the following samples:")
        print(donors_to_drop)
    df = pd.DataFrame(columns=[*adata_cell_pop.var_names, *obs_to_keep])

    adata_cell_pop.obs[donor_key] = adata_cell_pop.obs[donor_key].astype("category")
    
    for i, donor in enumerate(donors := adata_cell_pop.obs[donor_key].cat.categories):
        print(f"\tProcessing donor {i+1} out of {len(donors)}...", end="\r")
        
        if donor not in donors_to_drop:
            adata_donor = adata_cell_pop[adata_cell_pop.obs[donor_key] == donor]
            
            # create replicates for each donor
            indices = list(adata_donor.obs_names)
            random.shuffle(indices)
            indices = np.array_split(np.array(indices), replicates_per_patient)
            for i, rep_idx in enumerate(indices):
                adata_replicate = adata_donor[rep_idx]
                
                # specify how to aggregate: sum gene expression 
                # for each gene for each donor and also keep the condition information
                agg_dict = {gene: "sum" for gene in adata_replicate.var_names}
                for obs in obs_to_keep:
                    agg_dict[obs] = "first"
                    
                # create a df with all genes, donor and condition info
                df_donor = pd.DataFrame(adata_replicate.X.A)
                df_donor.index = adata_replicate.obs_names
                df_donor.columns = adata_replicate.var_names
                df_donor = df_donor.join(adata_replicate.obs[obs_to_keep])
                
                # aggregate
                df_donor = df_donor.groupby(donor_key).agg(agg_dict)
                df_donor[donor_key] = donor
                df.loc[f"{donor}_{i}"] = df_donor.loc[donor]
    print("\n")
    
    # create AnnData object from the df
    adata_cell_pop = sc.AnnData(
        df[adata_cell_pop.var_names], obs=df.drop(columns=adata_cell_pop.var_names)
    )
    return adata_cell_pop
