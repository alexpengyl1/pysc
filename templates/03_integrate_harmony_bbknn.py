# 03 Multi-sample integration: Harmony (preferred) / BBKNN / ingest
#
# [Purpose]
#   Removes batch effects between samples so that the same cell type clusters together whichever sample it comes from, then re-clusters and redraws the UMAP.
#   Three methods are provided; run only one:
#     A Harmony: corrects batches in PCA space (most common)
#     B BBKNN: builds a batch-balanced neighbour graph without changing the PCA
#     C ingest: when an annotated reference exists, projects the reference labels and UMAP onto the new data
# [Pipeline position] Upstream 02 -> this step -> downstream 05_annotate (alternative to 04 scVI)
# [Environment] sc-core (course ps7)
# [Input]
#   - results/02_qc_cluster/objects/02_processed.h5ad (needs obsm['X_pca'], obs['sample'])
#   - Method C also needs an annotated reference h5ad
# [Output] All in results/03_integrate/
#   - objects/03_integrated.h5ad: obsm['X_pca_harmony'] (method A), recomputed leiden and X_umap
#   - figures/: 03_umap_before_integration, 03_umap_after_harmony (or bbknn / ingest), 03_markers_dotplot_leiden
#   - tables/: 03_markers_leiden.csv, 03_cells_per_cluster_sample.csv (cells from each sample in each cluster)
# [How to read the results]
#   - UMAP before vs after: after integration, colouring by sample should look mixed while colouring by markers still separates cell types
#   - The same cell type still split by sample = under-correction; different cell types merged together = over-correction
#   - cells_per_cluster_sample: a cluster coming from only one sample is either residual batch effect or a genuine sample-specific population (e.g. a state seen only after treatment)

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/03_integrate/)
import pandas as pd
import scanpy as sc
import fair_utils as fair

IN_FILE = "results/02_qc_cluster/objects/02_processed.h5ad"
METHOD = "harmony"       # "harmony" / "bbknn" / "ingest"; for the record only, you choose which block to run
BATCH_KEY = "sample"
N_PCS = 30
RESOLUTION = 0.8
dirs = fair.setup_step("03_integrate")

# %%
# Step: read the previous result, plot the UMAP before integration
# Input: IN_FILE
# Output: adata; figures/03_umap_before_integration.pdf/.png
adata = sc.read_h5ad(IN_FILE)
sc.pl.umap(adata, color=[BATCH_KEY, "leiden"], show=False)
fair.save_fig(dirs, "umap_before_integration")

# %%
# Method A - Harmony
# Step: correct batches in PCA space
# Input: adata.obsm['X_pca']
# Output: adata.obsm['X_pca_harmony']
sc.external.pp.harmony_integrate(adata, key=BATCH_KEY, basis="X_pca", adjusted_basis="X_pca_harmony")

# %%
# Step: recompute neighbours, clusters and UMAP on the corrected space
# Input: adata.obsm['X_pca_harmony']
# Output: adata.obs['leiden'], adata.obsm['X_umap']; figures/03_umap_after_harmony.pdf/.png
sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_pcs=N_PCS, random_state=0)
sc.tl.leiden(adata, resolution=RESOLUTION, flavor="igraph", n_iterations=2, directed=False, random_state=0)
sc.tl.umap(adata, random_state=0)
sc.pl.umap(adata, color=[BATCH_KEY, "leiden"], show=False)
fair.save_fig(dirs, "umap_after_harmony")

# %%
# Method B - BBKNN (alternative to A: builds a batch-balanced neighbour graph directly, leaves the expression matrix and PCA unchanged)
# Step: BBKNN neighbour graph + clustering + UMAP
# Input: adata.obsm['X_pca']
# Output: adata.obs['leiden'], adata.obsm['X_umap']; figures/03_umap_after_bbknn.pdf/.png
# sc.external.pp.bbknn(adata, batch_key=BATCH_KEY, n_pcs=N_PCS)
# sc.tl.leiden(adata, resolution=RESOLUTION, flavor="igraph", n_iterations=2, directed=False, random_state=0)
# sc.tl.umap(adata, random_state=0)
# sc.pl.umap(adata, color=[BATCH_KEY, "leiden"], show=False)
# fair.save_fig(dirs, "umap_after_bbknn")

# %%
# Method C - ingest (with an annotated reference, projects the reference labels and UMAP onto the new data)
# Step: keep genes shared by both datasets, compute PCA/neighbours/UMAP on the reference, then map
# Input: adata_ref (annotated, obs['cell_type']), adata
# Output: adata.obs['cell_type'] (from the reference), adata.obsm['X_umap']; figures/03_umap_ingest_cell_type.pdf/.png
# adata_ref = sc.read_h5ad("data/reference.h5ad")
# genes = adata_ref.var_names.intersection(adata.var_names)
# adata_ref = adata_ref[:, genes].copy()
# adata = adata[:, genes].copy()
# sc.pp.pca(adata_ref)
# sc.pp.neighbors(adata_ref)
# sc.tl.umap(adata_ref)
# sc.tl.ingest(adata, adata_ref, obs="cell_type")
# sc.pl.umap(adata, color="cell_type", show=False)
# fair.save_fig(dirs, "umap_ingest_cell_type")

# %%
# Step: cells from each sample in each cluster (check for residual batch effects)
# Input: adata.obs
# Output: tables/03_cells_per_cluster_sample.csv
cluster_sample = pd.crosstab(adata.obs["leiden"], adata.obs[BATCH_KEY])
cluster_sample.to_csv(fair.out(dirs, "tables", "cells_per_cluster_sample.csv"))
cluster_sample

# %%
# Step: marker genes after integration
# Input: adata
# Output: tables/03_markers_leiden.csv; figures/03_markers_dotplot_leiden.pdf/.png
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=False)
sc.get.rank_genes_groups_df(adata, group=None).to_csv(fair.out(dirs, "tables", "markers_leiden.csv"), index=False)
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, dendrogram=False, show=False)
fair.save_fig(dirs, "markers_dotplot_leiden")

# %%
# Step: save
# Input: adata
# Output: objects/03_integrated.h5ad
adata.write(fair.out(dirs, "objects", "integrated.h5ad"))

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
