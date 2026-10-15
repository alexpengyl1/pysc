# 02 QC, normalisation, highly variable genes, PCA, clustering, UMAP, marker genes
#
# [Purpose]
#   Standard scanpy preprocessing of the merged raw data: remove low-quality cells and doublets, normalise, reduce dimensions and cluster,
#   and find marker genes for each cluster, ready for integration and annotation.
#   Workflow: QC metrics -> choose thresholds from plots -> filter -> Scrublet doublets -> normalise + log1p -> HVGs -> PCA -> neighbours -> Leiden -> UMAP -> markers
# [Pipeline position] Upstream 01 -> this step -> downstream 03 (Harmony/BBKNN) or 04 (scVI)
# [Environment] sc-core (course ps2 standard scanpy workflow, ps7)
# [Input]
#   - results/01_read_data/objects/01_merged_raw.h5ad (X = raw counts, obs['sample'])
# [Output] All in results/02_qc_cluster/
#   - objects/02_processed.h5ad: X = log1p-normalised (all genes, not scaled); layers['counts'] = raw counts; layers['log1p'] = same as X;
#     obs has QC metrics, predicted_doublet, leiden; obsm has X_pca, X_umap
#   - figures/: 02_qc_violin_by_sample, 02_qc_scatter_counts_genes, 02_hvg, 02_pca_variance_ratio, 02_umap_leiden_sample, 02_markers_dotplot_leiden (each as .pdf and .png)
#   - tables/: 02_qc_cells_kept.csv (cells per sample before/after filtering), 02_doublets_per_sample.csv, 02_markers_leiden.csv
# [How to read the results]
#   - QC violin plots: very low n_genes_by_counts = empty droplets or debris, very high = possible doublets; high pct_counts_mt = damaged cells. Set the thresholds at the top after looking at the plots
#   - qc_cells_kept table: a sample losing far more cells than the others has poor quality and needs a closer look
#   - The doublet rate usually rises with the number of cells loaded, roughly 1-8%
#   - PCA variance plot: PCs after the elbow add little; set N_PCS near the elbow
#   - UMAP coloured by sample: if samples form separate islands there is a batch effect, so integrate in the next step
#   - Marker dotplot: dot size = fraction of cells in the cluster expressing the gene, colour = mean expression; use it to decide what each leiden cluster is
# [Notes]
#   - Thresholds (MIN_GENES, MAX_PCT_MT etc.) depend on the tissue; do not copy the defaults blindly. Mouse mitochondrial genes start with mt-

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/02_qc_cluster/)
import pandas as pd
import scanpy as sc
import fair_utils as fair

IN_FILE = "results/01_read_data/objects/01_merged_raw.h5ad"
MT_PREFIX = "MT-"        # human "MT-", mouse "mt-"
MIN_GENES = 200          # course PBMC: 200-2500 genes, mito < 5%; tumour tissue: > 500 genes, mito < 20%
MAX_GENES = 6000
MAX_PCT_MT = 20
N_HVG = 2000
N_PCS = 30
RESOLUTION = 0.8
dirs = fair.setup_step("02_qc_cluster")

# %%
# Step: read the merged raw data
# Input: IN_FILE
# Output: adata
adata = sc.read_h5ad(IN_FILE)
adata

# %%
# Step: compute QC metrics (mitochondrial and ribosomal gene fractions)
# Input: adata
# Output: adata.obs gains n_genes_by_counts, total_counts, pct_counts_mt, pct_counts_ribo
adata.var["mt"] = adata.var_names.str.startswith(MT_PREFIX)
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL", "Rps", "Rpl"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)

# %%
# Step: QC violin plots by sample, to choose thresholds
# Input: adata.obs
# Output: figures/02_qc_violin_by_sample.pdf/.png
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], groupby="sample",
             multi_panel=True, show=False)
fair.save_fig(dirs, "qc_violin_by_sample")

# %%
# Step: total counts vs number of genes, coloured by mitochondrial fraction
# Input: adata.obs
# Output: figures/02_qc_scatter_counts_genes.pdf/.png
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt", show=False)
fair.save_fig(dirs, "qc_scatter_counts_genes")

# %%
# Step: filter cells and genes by the thresholds
# Input: adata
# Output: adata (filtered)
cells_before = adata.obs["sample"].value_counts()
sc.pp.filter_cells(adata, min_genes=MIN_GENES)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[(adata.obs.n_genes_by_counts < MAX_GENES) & (adata.obs.pct_counts_mt < MAX_PCT_MT)].copy()
adata

# %%
# Step: record cells per sample before and after filtering
# Input: cells_before, adata
# Output: tables/02_qc_cells_kept.csv
qc_kept = pd.DataFrame({"before": cells_before, "after": adata.obs["sample"].value_counts()})
qc_kept.to_csv(fair.out(dirs, "tables", "qc_cells_kept.csv"))
qc_kept

# %%
# Step: predict doublets per sample (Scrublet) and remove them
# Input: adata (raw counts)
# Output: adata (cells with predicted_doublet == True removed); tables/02_doublets_per_sample.csv
sc.pp.scrublet(adata, batch_key="sample", random_state=0)
doublets = pd.crosstab(adata.obs["sample"], adata.obs["predicted_doublet"])
doublets.to_csv(fair.out(dirs, "tables", "doublets_per_sample.csv"))
adata = adata[~adata.obs["predicted_doublet"]].copy()
doublets

# %%
# Step: keep the raw counts, then normalise total counts and log1p
# Input: adata (raw counts)
# Output: adata.X = log1p-normalised; layers['counts'], layers['log1p']
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log1p"] = adata.X.copy()

# %%
# Step: select highly variable genes from raw counts (seurat_v3, per sample)
# Input: adata.layers['counts']
# Output: adata.var['highly_variable']; figures/02_hvg.pdf/.png
sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat_v3", layer="counts", batch_key="sample")
sc.pl.highly_variable_genes(adata, show=False)
fair.save_fig(dirs, "hvg")

# %%
# Step: scale and PCA on the highly variable genes only (adata.X is not changed)
# Input: adata
# Output: adata.obsm['X_pca']
hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(hvg, max_value=10)
sc.pp.pca(hvg, n_comps=50, random_state=0)
adata.obsm["X_pca"] = hvg.obsm["X_pca"]
adata.uns["pca"] = hvg.uns["pca"]

# %%
# Step: PCA variance ratio plot, to choose N_PCS
# Input: adata.uns['pca']
# Output: figures/02_pca_variance_ratio.pdf/.png
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
fair.save_fig(dirs, "pca_variance_ratio")

# %%
# Step: neighbour graph, Leiden clustering, UMAP
# Input: adata.obsm['X_pca']
# Output: adata.obs['leiden'], adata.obsm['X_umap']; figures/02_umap_leiden_sample.pdf/.png
sc.pp.neighbors(adata, n_pcs=N_PCS, random_state=0)
sc.tl.leiden(adata, resolution=RESOLUTION, flavor="igraph", n_iterations=2, directed=False, random_state=0)
sc.tl.umap(adata, random_state=0)
sc.pl.umap(adata, color=["leiden", "sample", "pct_counts_mt"], show=False)
fair.save_fig(dirs, "umap_leiden_sample")

# %%
# Step: marker genes for each leiden cluster (Wilcoxon)
# Input: adata (X = log1p)
# Output: tables/02_markers_leiden.csv; figures/02_markers_dotplot_leiden.pdf/.png
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=False)
sc.get.rank_genes_groups_df(adata, group=None).to_csv(fair.out(dirs, "tables", "markers_leiden.csv"), index=False)
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, dendrogram=False, show=False)
fair.save_fig(dirs, "markers_dotplot_leiden")

# %%
# Step: save the processed data
# Input: adata
# Output: objects/02_processed.h5ad
adata.write(fair.out(dirs, "objects", "processed.h5ad"))

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
