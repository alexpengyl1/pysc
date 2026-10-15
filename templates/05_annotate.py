# 05 Cell type annotation: manual marker annotation + CellTypist cross-check
#
# [Purpose]
#   Names the cell type of each cluster. First look at known marker genes across clusters, then annotate automatically with a pretrained CellTypist model,
#   compare the two and write a cluster -> cell type table by hand. This creates obs['cell_type'], which every later analysis uses.
#   Workflow: marker dotplot -> CellTypist annotation -> compare with clusters -> fill in the mapping by hand -> plot annotation -> save
# [Pipeline position] Upstream 03 or 04 -> this step -> downstream: every analysis in 06-24 (all read 05_annotated.h5ad)
# [Environment] sc-core (course ps2/ps6 manual annotation, ps9 CellTypist)
# [Input]
#   - results/03_integrate/objects/03_integrated.h5ad (or the output of 02 / 04); X = log1p (1e4 per cell, all genes)
#   - MARKERS at the top: known markers for the cell types in your tissue (genes missing from the data are dropped and printed)
#   - CellTypist model: downloaded on first run (needs internet)
# [Output] All in results/05_annotate/
#   - objects/05_annotated.h5ad: obs['cell_type'] (final annotation), obs['celltypist_*'] (automatic labels and confidence)
#   - figures/: 05_markers_dotplot_leiden, 05_umap_celltypist, 05_celltypist_vs_leiden_dotplot, 05_umap_cell_type, 05_markers_dotplot_cell_type
#   - tables/: 05_cluster_annotation.csv (cluster -> cell type), 05_celltypist_vs_leiden.csv, 05_cells_per_type_sample.csv
# [How to read the results]
#   - Marker dotplot: a set of markers all high in one cluster (large, dark dots) = that cell type
#   - CellTypist: majority_voting (per-cluster vote) is more stable than per-cell predicted_labels; recheck cells with low conf_score
#   - celltypist_vs_leiden: one cluster matching several CellTypist labels may be under-clustered; where the two disagree, go back to the markers
#   - cells_per_type_sample: cell composition per sample, the raw data for the composition analysis in 06
# [Notes]
#   - Fill in the mapping table after looking at the plots above; do not use the defaults as they are
#   - Most CellTypist models are human; Adult_Mouse_Gut works for mouse gut, be careful with other mouse data

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/05_annotate/)
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models
import fair_utils as fair

IN_FILE = "results/03_integrate/objects/03_integrated.h5ad"
CLUSTER_KEY = "leiden"
CELLTYPIST_MODEL = "Immune_All_Low.pkl"    # others: Immune_All_High, Cells_Intestinal_Tract, Adult_Mouse_Gut ...
MARKERS = {                                # adapt to your tissue
    "T": ["CD3D", "CD3E"],
    "B": ["MS4A1", "CD79A"],
    "Myeloid": ["CD14", "LYZ"],
    "Epithelial": ["EPCAM", "KRT8"],
    "Fibroblast": ["COL1A1", "DCN"],
    "Endothelial": ["PECAM1", "VWF"],
}
dirs = fair.setup_step("05_annotate")

# %%
# Step: read the data, keep only marker genes present in the data
# Input: IN_FILE, MARKERS
# Output: adata, markers (dict without missing genes); prints the dropped genes
adata = sc.read_h5ad(IN_FILE)
markers = {cell: [g for g in genes if g in adata.var_names] for cell, genes in MARKERS.items()}
markers = {cell: genes for cell, genes in markers.items() if genes}
print("Markers not found in the data:", [g for genes in MARKERS.values() for g in genes if g not in adata.var_names])

# %%
# Step: marker dotplot by leiden cluster
# Input: adata, markers
# Output: figures/05_markers_dotplot_leiden.pdf/.png
sc.pl.dotplot(adata, markers, groupby=CLUSTER_KEY, show=False)
fair.save_fig(dirs, "markers_dotplot_leiden")

# %%
# Step: download CellTypist models (needed on first run), list available models
# Input: none
# Output: local model cache; printed model descriptions
models.download_models(force_update=False)
models.models_description()

# %%
# Step: CellTypist annotation (majority voting)
# Input: adata (X = log1p)
# Output: predictions; adata.obs['celltypist_label'], ['celltypist_majority'], ['celltypist_conf']
predictions = celltypist.annotate(adata, model=CELLTYPIST_MODEL, majority_voting=True)
adata.obs["celltypist_label"] = predictions.predicted_labels["predicted_labels"]
adata.obs["celltypist_majority"] = predictions.predicted_labels["majority_voting"]
adata.obs["celltypist_conf"] = predictions.probability_matrix.max(axis=1)

# %%
# Step: CellTypist results on the UMAP
# Input: adata.obs
# Output: figures/05_umap_celltypist.pdf/.png
sc.pl.umap(adata, color=[CLUSTER_KEY, "celltypist_majority", "celltypist_conf"], legend_loc="on data", show=False)
fair.save_fig(dirs, "umap_celltypist")

# %%
# Step: CellTypist results vs clusters (plot + cross table)
# Input: predictions, adata.obs
# Output: figures/05_celltypist_vs_leiden_dotplot.pdf/.png; tables/05_celltypist_vs_leiden.csv
celltypist.dotplot(predictions, use_as_reference=CLUSTER_KEY, use_as_prediction="majority_voting", show=False)
fair.save_fig(dirs, "celltypist_vs_leiden_dotplot")
pd.crosstab(adata.obs[CLUSTER_KEY], adata.obs["celltypist_majority"]).to_csv(fair.out(dirs, "tables", "celltypist_vs_leiden.csv"))

# %%
# Step: write the cluster -> cell type mapping by hand (after looking at the plots above) and save it
# Input: adata.obs[CLUSTER_KEY]
# Output: adata.obs['cell_type']; tables/05_cluster_annotation.csv
annotation = {
    "0": "T",
    "1": "B",
    # ...
}
adata.obs["cell_type"] = adata.obs[CLUSTER_KEY].map(annotation).astype("category")
pd.Series(annotation, name="cell_type").rename_axis(CLUSTER_KEY).to_csv(fair.out(dirs, "tables", "cluster_annotation.csv"))

# %%
# Step: UMAP and marker dotplot of the annotation
# Input: adata
# Output: figures/05_umap_cell_type.pdf/.png, 05_markers_dotplot_cell_type.pdf/.png
sc.pl.umap(adata, color="cell_type", legend_loc="on data", show=False)
fair.save_fig(dirs, "umap_cell_type")
sc.pl.dotplot(adata, markers, groupby="cell_type", show=False)
fair.save_fig(dirs, "markers_dotplot_cell_type")

# %%
# Step: number of cells of each cell type per sample
# Input: adata.obs
# Output: tables/05_cells_per_type_sample.csv
pd.crosstab(adata.obs["cell_type"], adata.obs["sample"]).to_csv(fair.out(dirs, "tables", "cells_per_type_sample.csv"))

# %%
# Step: save
# Input: adata
# Output: objects/05_annotated.h5ad
adata.write(fair.out(dirs, "objects", "annotated.h5ad"))

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
