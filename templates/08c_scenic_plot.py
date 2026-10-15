# 08c Find cell-type-specific regulons and plot them
#
# [Purpose]
#   Reads the regulon activity matrix from 08b back into AnnData, finds the transcription factors specifically active in each cell type, and draws a matrix plot and a clustered heatmap.
# [Pipeline position] Upstream 08b -> this step (end point)
# [Environment] sc-core (course ps13; the course only drew the clustermap, a per-cell-type summary is added here)
# [Input]
#   - results/08_scenic/tables/08_auc_matrix.csv (from 08b)
#   - results/08_scenic/tables/08_cell_info.csv (from 08a)
# [Output] All in results/08_scenic/
#   - tables/08_regulon_markers_cell_type.csv: regulons specific to each cell type (Wilcoxon, one vs rest)
#   - figures/: 08_regulon_matrixplot_cell_type, 08_auc_clustermap
# [How to read the results]
#   - Matrix plot: one row per cell type, one column per regulon, colour = mean activity (scaled per column); a bright regulon in one row = an identity or state regulator of that cell type
#   - regulon_markers table: regulons with high scores and small pvals_adj are the most specific
#   - Clustermap: cells clustered by regulon activity; whether they separate by cell type reflects differences in regulatory programmes

# %%
# Step: import packages, set paths and output folder
# Input: none
# Output: dirs (results/08_scenic/)
import pandas as pd
import seaborn as sns
import scanpy as sc
import anndata as ad
import fair_utils as fair

AUC_FILE = "results/08_scenic/tables/08_auc_matrix.csv"
CELL_INFO_FILE = "results/08_scenic/tables/08_cell_info.csv"
N_TOP = 3
dirs = fair.setup_step("08_scenic")

# %%
# Step: read the AUC matrix into an AnnData with regulons as features
# Input: AUC_FILE, CELL_INFO_FILE
# Output: auc_adata (cells x regulons, obs['cell_type'])
auc = pd.read_csv(AUC_FILE, index_col=0)
cell_info = pd.read_csv(CELL_INFO_FILE, index_col=0)
auc_adata = ad.AnnData(auc)
auc_adata.obs["cell_type"] = cell_info.loc[auc.index, "cell_type"].astype("category")
auc_adata

# %%
# Step: regulons specific to each cell type (Wilcoxon, one vs rest)
# Input: auc_adata
# Output: tables/08_regulon_markers_cell_type.csv
sc.tl.rank_genes_groups(auc_adata, "cell_type", method="wilcoxon")
sc.get.rank_genes_groups_df(auc_adata, group=None).to_csv(
    fair.out(dirs, "tables", "regulon_markers_cell_type.csv"), index=False)

# %%
# Step: matrix plot of the top N_TOP specific regulons per cell type
# Input: auc_adata
# Output: figures/08_regulon_matrixplot_cell_type.pdf/.png
sc.pl.rank_genes_groups_matrixplot(auc_adata, n_genes=N_TOP, standard_scale="var", dendrogram=False, show=False)
fair.save_fig(dirs, "regulon_matrixplot_cell_type")

# %%
# Step: clustered heatmap of all cells x regulons (the plot from course ps13)
# Input: auc
# Output: figures/08_auc_clustermap.pdf/.png
grid = sns.clustermap(auc, figsize=(12, 12), yticklabels=False)
fair.save_fig(dirs, "auc_clustermap", fig=grid)

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info_08c.json, logs/pip_freeze_08c.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()}, script="08c")
