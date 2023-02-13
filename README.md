# Notebooks and scripts for epicardioids analysis
* `notebooks`
  * `analyse_signac_unfiltered.ipynb`: QC filtering of ATAC data
  * `cellrank_analysis.ipynb`: CellRank analysis on metacells
  * `create_anndata.ipynb`: Loading of raw data into anndata
  * `epicardial_subclustering.ipynb`: Analysis for epicardial cells
  * `epicardioids_metacell_analysis.ipynb `: Analysis of metacells, clustering, DEG
  * `impute_weights_magic.ipynb `: Gene activity imputation on metacell embedding
  * `jcf_subclustering.ipynb`: Analysis for JCF cells
  * `pando_grn.ipynb `: Pando GRN inference on JCF cells
  * `scglue_unfiltered.ipynb `: scGLUE to generate metacell embeddings (this needs bipartite matching)
  * `scvelo_analysis.ipynb`: scvelo analysis on RNA modality

* `scripts`
  * `amulet`: Doublet detection using Amulet
  * `integration`: Bipartite matching (matching ATAC to RNA cells) on scGLUE embedding 
