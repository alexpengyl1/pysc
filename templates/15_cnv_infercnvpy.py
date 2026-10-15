# 15 inferCNVpy: infer copy number variation, separate malignant from normal cells
#
# [Purpose]
#   Smooths the expression of neighbouring genes along each chromosome and compares it with reference (normal) cells
#   to infer which chromosome segments are gained or lost in each cell. Large CNVs are typical of tumour cells and help separate malignant from normal cells.
#   Workflow: add gene positions from GTF -> infer CNV against normal cells -> chromosome heatmap -> cluster and score in CNV space -> mark tumour clusters -> save
# [Pipeline position] Upstream 05 -> this step (end point)
# [Environment] sc-core (course ps20)
# [Input]
#   - results/05_annotate/objects/05_annotated.h5ad (X = log1p; obs['cell_type'])
#   - A GTF matching the reference genome used for alignment (gives genes their chromosome positions)
#   - REFERENCE_TYPES at the top: cell types known to be normal (e.g. immune cells, fibroblasts)
# [Output] All in results/15_infercnv/
#   - objects/15_infercnv.h5ad: obsm['X_cnv'], obs['cnv_score'], ['cnv_leiden'], ['cnv_status']
#   - tables/: 15_cnv_score_per_cell.csv, 15_cnv_status_by_cell_type.csv
#   - figures/: 15_chromosome_heatmap_cell_type, 15_chromosome_heatmap_cnv_leiden, 15_cnv_umap, 15_umap_cnv_status
# [How to read the results]
#   - Chromosome heatmap: reference cells should be mostly flat; red = gain, blue = loss
#   - Large segmental gains / losses concentrated in epithelial or tumour clusters -> malignant; cnv_leiden clusters with high cnv_score
#   - cnv_status_by_cell_type: number of cells called tumour in each cell type
# [Notes]
#   - Fill in TUMOR_CLUSTERS by hand after looking at the heatmap
#   - Non-tumour samples (e.g. radiation-injured tissue) usually have no large CNVs; the result only helps rule out malignant transformation

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/15_infercnv/)
import pandas as pd
import scanpy as sc
import infercnvpy as cnv
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
GTF = "data/genes.gtf"
REFERENCE_TYPES = ["T", "B", "Myeloid", "Endothelial"]   # cell types used as the normal reference
WINDOW_SIZE = 250
TUMOR_CLUSTERS = ["1", "4"]                               # fill in after looking at the cnv_leiden heatmap
dirs = fair.setup_step("15_infercnv")

# %%
# Step: read the data, add genomic positions from the GTF
# Input: IN_FILE, GTF
# Output: adata.var gains chromosome, start, end
adata = sc.read_h5ad(IN_FILE)
cnv.io.genomic_position_from_gtf(GTF, adata, gtf_gene_id="gene_name")
adata.var[["chromosome", "start", "end"]].head()

# %%
# Step: infer CNV against the normal reference cells
# Input: adata
# Output: adata.obsm['X_cnv']
cnv.tl.infercnv(adata, reference_key="cell_type", reference_cat=REFERENCE_TYPES, window_size=WINDOW_SIZE)

# %%
# Step: chromosome heatmap by cell type
# Input: adata.obsm['X_cnv']
# Output: figures/15_chromosome_heatmap_cell_type.pdf/.png
cnv.pl.chromosome_heatmap(adata, groupby="cell_type", show=False)
fair.save_fig(dirs, "chromosome_heatmap_cell_type")

# %%
# Step: cluster in CNV space and compute a CNV score per cell
# Input: adata.obsm['X_cnv']
# Output: obs['cnv_leiden'], obs['cnv_score'], obsm['X_cnv_umap']
cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata)
cnv.tl.umap(adata)
cnv.tl.cnv_score(adata)

# %%
# Step: chromosome heatmap by CNV cluster
# Input: adata
# Output: figures/15_chromosome_heatmap_cnv_leiden.pdf/.png
cnv.pl.chromosome_heatmap(adata, groupby="cnv_leiden", dendrogram=True, show=False)
fair.save_fig(dirs, "chromosome_heatmap_cnv_leiden")

# %%
# Step: UMAP in CNV space (CNV clusters, score, cell type)
# Input: adata
# Output: figures/15_cnv_umap.pdf/.png
cnv.pl.umap(adata, color=["cnv_leiden", "cnv_score", "cell_type"], show=False)
fair.save_fig(dirs, "cnv_umap")

# %%
# Step: mark clusters with clear CNVs as tumour
# Input: adata.obs['cnv_leiden'], TUMOR_CLUSTERS
# Output: obs['cnv_status']; figures/15_umap_cnv_status.pdf/.png
adata.obs["cnv_status"] = "normal"
adata.obs.loc[adata.obs["cnv_leiden"].isin(TUMOR_CLUSTERS), "cnv_status"] = "tumor"
sc.pl.umap(adata, color=["cnv_status", "cell_type"], show=False)
fair.save_fig(dirs, "umap_cnv_status")

# %%
# Step: save per-cell CNV results and tumour/normal counts per cell type
# Input: adata.obs
# Output: tables/15_cnv_score_per_cell.csv, 15_cnv_status_by_cell_type.csv
adata.obs[["cell_type", "cnv_leiden", "cnv_score", "cnv_status"]].to_csv(fair.out(dirs, "tables", "cnv_score_per_cell.csv"))
pd.crosstab(adata.obs["cell_type"], adata.obs["cnv_status"]).to_csv(fair.out(dirs, "tables", "cnv_status_by_cell_type.csv"))

# %%
# Step: save
# Input: adata
# Output: objects/15_infercnv.h5ad
adata.write(fair.out(dirs, "objects", "infercnv.h5ad"))

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
