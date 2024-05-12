import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns
from anndata import AnnData
import os


def plot_qc_distributions(adata: AnnData, save_fig_to: str = None,
                          suffix: str = '') -> None:
    """Plot quality control distributions of UMI counts, expressed genes,
    and mitochondrial gene expression.
    Parameters:
        adata (AnnData):
            Anndata object containing the data.
        save_fig_to (str):
            Directory to save the figure. None by default (doesn't save to file).
        suffix (str):
            ('before', 'after')
            If `save_fig_to` passed, suffix to add at the end of the
            filename before extension
    """

    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(15, 3), dpi=100)

    # Draw the distribution of the total number of UMIs
    # take a logarithm to (visually) reduce the spread of data
    sns.histplot(np.log(adata.obs["total_counts"]), kde=False, ax=axs[0])
    axs[0].set_xlabel("log(UMI per cell)")
    axs[0].set_ylabel("Number of cells")
    axs[0].set_title("UMI distribution")

    # Draw the distribution of the n genes per cell
    sns.histplot(adata.obs["n_genes_by_counts"], kde=False, ax=axs[1])
    axs[1].set_xlabel("Number of expressed genes per cell")
    axs[1].set_ylabel("Number of cells")
    axs[1].set_title("Distribution of expressed genes")

    # Draw the distribution of the mitochondrial expression
    sns.histplot(adata.obs["pct_counts_mt"], kde=False, ax=axs[2])
    axs[2].set_xlabel("% of mitochondrial genes per cell")
    axs[2].set_ylabel("Number of cells")
    axs[2].set_title("Distribution of mt genes")

    fig.tight_layout()

    if save_fig_to:
        plt.savefig(
            os.path.join(save_fig_to, f'qc_{suffix}'),
            dpi=300,
            bbox_inches='tight'
        )


def plot_scrublet_distribution(adata: AnnData, save_fig_to: str = None) -> None:
    """
    Plot the distribution of doublet scores and simulated doublet scores
    from Scrublet analysis.

    Parameters:
    adata (AnnData):
        Annotated data object containing doublet scores and simulated doublet scores.
    save_fig_to (str):
        Directory to save the figure. None by default (doesn't save to file).

    Returns:
    None
    """
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.kdeplot(adata.obs["doublet_score"], color="blue", ax=ax)
    sns.kdeplot(adata.uns["scrublet"]["doublet_scores_sim"], color="red", ax=ax)
    ax.plot([adata.uns["scrublet"]["threshold"], adata.uns["scrublet"]["threshold"]], [0, 20])
    ax.set_ylim(0, 20)
    ax.grid(alpha=0.3)

    if save_fig_to:
        plt.savefig(
            os.path.join(save_fig_to, 'scrublet'),
            dpi=300,
            bbox_inches='tight'
        )


def barplot(adata: AnnData, cluster_key: str, batch_key: str,
            save_fig_to: str = None) -> None:
    """
    Creates a barplot with the distribution of cells
    from different batches across clusters.

    Parameters:
    adata : AnnData
        Annotated data matrix.
    cluster_key : str
        Name of column in `adata.obs` that corresponds to the name of clusters.
    batch_key : str
        Name of column in `adata.obs`
        that corresponds to the name of samples / batches.
    save_fig_to: str
        Directory to save the figure. None by default (doesn't save to file).

    Returns:
    None.
    """
    fig_width = len(adata.obs[cluster_key].cat.categories) * 0.3
    fig, ax = plt.subplots(dpi=150, figsize=(fig_width, 2))

    sizes = adata.obs.groupby([batch_key, cluster_key]).size()
    props = sizes.groupby(level=1).apply(lambda x: 100 * x / x.sum()).reset_index(allow_duplicates=True,
                                                                                 # level=[1, 2]
                                                                                 )
    props = props.pivot(columns=cluster_key, index=batch_key).T
    props.index = props.index.droplevel(0)
    props.fillna(0, inplace=True)
    props.plot.bar(stacked=True, width=1,
                   edgecolor="black", ax=ax
                   )
    plt.xticks(rotation=90)
    ax.set_xlabel("")
    ax.legend(loc=(1.01, 0.45), edgecolor="white")

    if save_fig_to:
        plt.savefig(
            os.path.join(save_fig_to, f'{batch_key}_barplot'),
            dpi=300,
            bbox_inches='tight'
        )

def plot_leiden_clusters(adata: AnnData,
                         save_fig_to: str = None,
                         suffix: str = '',
                         leiden_key: str = 'leiden') -> None:
    """Plot UMAP and dendrogram of Leiden clusters
    Args:
        adata (AnnData):
            Annotated data object
        save_fig_to (str):
            Directory to save the figure (optional)
        suffix (str):
            (e.g. 'all', 'CD8', etc.)
            If `save_fig_to` passed, suffix to add at the end of the
            filename before extension
        leiden_key (str): 
            Key in adata.obs with Leiden clusters
    Returns:
        None
    """
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    sc.pl.umap(adata, color=leiden_key, frameon=False, show=False, ax=ax[1],
               legend_loc="on data", legend_fontoutline=2)
    ax[1].set_title('')
    sc.pl.dendrogram(adata, groupby=leiden_key, orientation="left", show=False, ax=ax[0])

    if save_fig_to:
        plt.savefig(
            os.path.join(save_fig_to, f'leiden_{suffix}.png'),
            dpi=300,
            bbox_inches='tight'
        )
