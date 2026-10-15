# 19 Hotspot: find locally co-varying gene modules within one cell type
#
# [Purpose]
#   Within one cell type, finds genes whose expression changes consistently among similar cells, and groups them into modules.
#   Each module usually corresponds to a subtype or cell state (e.g. activation, stress, proliferation).
#   Workflow: subset cell type -> recompute PCA -> KNN graph -> gene autocorrelation -> local correlation of significant genes -> modules -> module scores -> UMAP
# [Pipeline position] Upstream 05 -> this step (end point; change CELLTYPE to repeat for each cell type)
# [Environment] sc-core (course ps24)
# [Input]
#   - One cell type from results/05_annotate/objects/05_annotated.h5ad (layers['counts'], obs['total_counts']; latent space = recomputed X_pca, X_scVI also works)
# [Output] All in results/19_hotspot/, file names include the cell type
#   - tables/: 19_autocorrelation_<celltype>.csv, 19_modules_<celltype>.csv, 19_module_scores_<celltype>.csv
#   - figures/: 19_local_correlations_<celltype>, 19_umap_module_scores_<celltype>
# [How to read the results]
#   - Autocorrelation table: high Z, small FDR = expression consistent among neighbouring cells, an informative gene for cell state
#   - Local correlation heatmap: blocks on the diagonal = modules; module -1 = not assigned to any module
#   - modules table: the module of each gene; run enrichment on module genes to name the modules
#   - Module score UMAP: regions where a module is active, corresponding to subtypes or activation states

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/19_hotspot/)
import scanpy as sc
import hotspot
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
CELLTYPE = "T"
N_PCS = 10
N_NEIGHBORS = 30
N_GENES = 500             # use the top N significant genes by autocorrelation for modules
MIN_MODULE_GENES = 15
FDR = 0.05
N_JOBS = 4
SEED = 0
TAG = CELLTYPE.lower().replace(" ", "_")
dirs = fair.setup_step("19_hotspot")

# %%
# Step: read the data, take one cell type, recompute PCA
# Input: IN_FILE
# Output: adata (subset, obsm['X_pca'], layers['counts_csc'])
adata = sc.read_h5ad(IN_FILE)
adata = adata[adata.obs["cell_type"] == CELLTYPE].copy()
sc.pp.pca(adata, n_comps=N_PCS, mask_var="highly_variable", random_state=SEED)
adata.layers["counts_csc"] = adata.layers["counts"].tocsc()

# %%
# Step: build the Hotspot object (raw counts, negative binomial model) and the KNN graph
# Input: adata
# Output: hs
hs = hotspot.Hotspot(adata, layer_key="counts_csc", model="danb", latent_obsm_key="X_pca",
                     umi_counts_obs_key="total_counts")
hs.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS)

# %%
# Step: autocorrelation of each gene
# Input: hs
# Output: autocorr; tables/19_autocorrelation_<celltype>.csv
autocorr = hs.compute_autocorrelations(jobs=N_JOBS)
autocorr.to_csv(fair.out(dirs, "tables", f"autocorrelation_{TAG}.csv"))
autocorr.head(15)

# %%
# Step: local correlations between significant genes, grouped into modules
# Input: autocorr
# Output: modules; tables/19_modules_<celltype>.csv
genes = autocorr.loc[autocorr.FDR < FDR].sort_values("Z", ascending=False).head(N_GENES).index
hs.compute_local_correlations(genes, jobs=N_JOBS)
modules = hs.create_modules(min_gene_threshold=MIN_MODULE_GENES, core_only=True, fdr_threshold=FDR)
modules.to_csv(fair.out(dirs, "tables", f"modules_{TAG}.csv"))
modules.value_counts()

# %%
# Step: local correlation heatmap
# Input: hs
# Output: figures/19_local_correlations_<celltype>.pdf/.png
hs.plot_local_correlations(vmin=-12, vmax=12)
fair.save_fig(dirs, f"local_correlations_{TAG}")

# %%
# Step: module scores
# Input: hs
# Output: module_scores; tables/19_module_scores_<celltype>.csv
module_scores = hs.calculate_module_scores()
module_scores.to_csv(fair.out(dirs, "tables", f"module_scores_{TAG}.csv"))

# %%
# Step: module scores on the UMAP
# Input: module_scores
# Output: figures/19_umap_module_scores_<celltype>.pdf/.png
for module in module_scores.columns:
    adata.obs[f"Module {module}"] = module_scores[module]
sc.pl.umap(adata, color=[f"Module {m}" for m in module_scores.columns], vmin=-1, vmax=1, show=False)
fair.save_fig(dirs, f"umap_module_scores_{TAG}")

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
