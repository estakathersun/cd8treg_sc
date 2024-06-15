# CD8<sup>+</sup>Treg scRNA-seq analysis
The repository contains code of scRNA-seq data analysis of the CD8+HLA-DR+ regulatory T cells as part of a research project 
**"Transcriptional signatures and age-related changes in CD8+HLA-DR+ regulatory T cells"** at the Bioinformatics Institute 2023/2024

## Introduction
First described in 2014, the CD8<sup>+</sup>HLA-DR<sup>+</sup> regulatory T lymphocytes (CD8<sup>+</sup>Treg) subset is known for its role in suppressing effector T cells through checkpoint inhibitor molecules, sharing certain features with conventional CD4<sup>+</sup>Treg cells. Despite this, the detailed nature and function of these cells are still not well understood. Understanding this subset is particularly significant in light of age-related alterations in the immune system and the heightened vulnerability of CD8<sup>+</sup> T lymphocytes to such changes. Preliminary findings indicated that the CD8<sup>+</sup>HLA-DR<sup>+</sup> population differentiates into two subpopulations based on CD127 (*IL7R*) surface expression. This led to defining the CD8<sup>+</sup>Treg phenotype as CD3<sup>+</sup>CD8<sup>+</sup>HLA-DR<sup>+</sup>CD127<sup>low</sup>. This research aimed to identify transcriptional signatures and examine age-related changes in gene expression within the CD8<sup>+</sup>Treg population, utilizing publicly accessible single-cell RNA-seq (scRNA-seq) data. 

## Workflow
**Data**: Single-cell RNA sequencing (scRNA-seq) data from peripheral blood mononuclear cells (PBMCs) of 982 donors, aged 19 to 97, from the [OneK1K cohort](https://doi.org/10.1126/science.abf3041) and 99 healthy donors, aged 22 to 75, from the [SLE study](https://doi.org/10.1126/science.abf1970) were analyzed.

The `scanpy` package was used for data analysis:
1. **Data preprocessing**: 
   - Normalization;
   - Logarithmization;
   - Scaling;
   - Dimensionality reduction (via PCA);
   - Batch correction (`harmony`);
   - Leiden clustering, UMAP representation.
2. **Automatic cell types annotation** with `celltypist` was employed to annotate the cell clusters and identify CD8<sup>+</sup> effector T cells.
3. **CD8<sup>+</sup> T cell subsets selection.**
4. **Return raw counts to the CD8<sup>+</sup> T cells** to better represent them in reduced dimensionality embedding.
5. **Repeat preprocessing steps for the CD8<sup>+</sup> T cells data.**
6. **Manual annotation of CD8<sup>+</sup>Treg cells** based on increased expression of *HLA-DRA*, *HLA-DRB1*, and *HLA-DRB5* genes and decreased expression of *IL7R* gene.
7. **Differential gene expression analysis (DEA) of CD8<sup>+</sup>Treg at the single-cell level.**
8. **DEA and gene set enrichment analysis (GSEA) of CD8<sup>+</sup>Treg at the pseudobulk level** with `edgeR`, `clusterProfiler` and `gprofiler2` **R** packages.

The workflow scheme is represented below:

<img width="800" src="https://github.com/estakathersun/cd8treg_sc/assets/143891764/6a52bbe0-3ea0-4db6-bbee-fd734bb27ef3">

## Results
###  1. The CD8<sup>+</sup>Treg transcriptional signatures
To determine CD8<sup>+</sup>Treg signatures, scRNA-seq data from donors aged 19 to 60 years (n = 299) were analyzed. According to previous experimental results, T cells in this age range typically do not exhibit distinct age-related changes. 

It should be noted that cells potentially corresponding to CD8<sup>+</sup>Treg **were not separated into a distinct cluster** using the Leiden algorithm but were distributed among CD8<sup>+</sup> effector T cells. Therefore, the CD8<sup>+</sup>Treg population was defined by their phenotype as CD8<sup>+</sup> effector T cells **with increased expression of *HLA-DRA*, *HLA-DRB1*, and *HLA-DRB5* genes and decreased expression of *IL7R* gene (CD127)**:

<img width="800" src="https://github.com/estakathersun/cd8treg_sc/assets/143891764/fc232221-1b7d-4005-8ff1-7c61187ddf28">

DEA at the single-cell level revealed that the CD8<sup>+</sup>Treg population, compared to the main population of CD8<sup>+</sup> effector T cells, exhibited **increased expression of genes associated with MHC-II-mediated antigen presentation, cytotoxicity, and cytoskeletal organization**. At the pseudobulk level, GSEA based on the GO terms also characterized the CD8<sup>+</sup>Treg population by increased expression of genes involved in antigen presentation, particularly those mediated by MHC class II molecules.

The single-cell-level DEA and pseudobulk-level GSEA results are represent below respectively:

<img width="500" align='top' src="https://github.com/estakathersun/cd8treg_sc/assets/143891764/c8e91f32-f3e8-4c97-8084-0bbf9591a506"><img width="500" align='top' src="https://github.com/estakathersun/cd8treg_sc/assets/143891764/1af99378-7c3e-4f2e-b56e-ed7095f96d25">
<br/><br/>

### 2. Age-related changes in the CD8<sup>+</sup>Treg subset
To assess age-associated changes at the transcriptomic level, CD8<sup>+</sup>Treg populations from young (20-35 years, n = 152) and old (70-97 years, n = 424) donors were analyzed.
DEA at the pseudobulk level among older donors showed a shift in the CD8<sup>+</sup>Treg transcriptional profile **towards a terminally differentiated phenotype**. Specifically, a decrease in the expression of cytotoxic molecules such as *LYZ*, *GZMA*, and *GZMK* was observed, while the expression of *GZMH*, characteristic of terminally differentiated CD8<sup>+</sup> T lymphocytes, increased:

<img width="300" src="https://github.com/estakathersun/cd8treg_sc/assets/143891764/1d36dbb3-b894-4afb-96f3-c935241b545e">

<br/><br/>
Additionally, CD8<sup>+</sup>Treg were characterized by increased expression levels of MHC molecules, as confirmed by GSEA:

<img width="500" align='top' src="https://github.com/estakathersun/cd8treg_sc/assets/143891764/7420d632-1769-4653-ac25-7776463321c2">

<br/><br/>
Furthermore, GSEA results indicated a **decrease** in the expression of genes **involved in RNA biosynthesis and metabolism, as well as T cell differentiation, activation**, and immune response with age. However, no increase in the expression of genes associated with exhaustion, such as PD-1, Tim3, Lag3, TIGIT, CD160, and CD244, was observed. It can be inferred that with age, the CD8<sup>+</sup>Treg population transitions to a terminally differentiated phenotype and exhibits a decreased functional response capacity, but this population does not undergo cellular exhaustion.

<img width="600" align='top' src="https://github.com/estakathersun/cd8treg_sc/assets/143891764/2c95c53a-13dc-4ffb-8f38-d38d0bbb161d">
<br/><br/>

## Conclusion
Using the scRNA-seq data, we observed that the CD8<sup>+</sup>Treg subpopulation is a heterogeneous group of CD8<sup>+</sup> effector T lymphocytes. This subpopulation shows increased expression of genes associated with cytotoxicity, cytoskeletal rearrangement, and MHC-II-mediated antigen presentation. This would suggest that CD8<sup>+</sup>Treg-mediated suppression likely involves cell-contact dependent cytolysis of target cells. In older adults, CD8<sup>+</sup>Treg cells exhibit changes in marker gene expression indicative of a terminally differentiated cell phenotype. Despite this shift, there is no evidence of increased expression of exhaustion markers.  However, there is a decrease in the expression of genes regulating key processes of T cell activation and function, suggesting that the suppressor function of CD8<sup>+</sup>Treg decreases with age.

## Contributors
K. Matveeva<sup>1,2</sup>, S. Kolmykov<sup>2</sup>, D. Shevyrev<sup>2</sup>

- <sup>1</sup> Bioinformatics Institute, Kantemirovskaya st. 2A, 197342, St. Petersburg, Russia
- <sup>2</sup> Sirius University of Science and Technology, Olympic Ave., 1, 354340, Sirius, Russia

## Contacts
For any questions or suggestions, please feel free to reach out to [Kseniia Matveeva](https://github.com/estakathersun) at [email](mailto:ksum.miha@gmail.com) or [Telegram](https://t.me/estakathersun)✨

