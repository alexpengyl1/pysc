# 24 CellOracle: build gene regulatory networks (GRNs) per cluster and compare network centrality
#
# [Purpose]
#   Combines a prior network of TF binding sites (base GRN) with single-cell expression to build a TF -> target network for each cluster,
#   then scores each gene's importance in the network (centrality) and compares key regulators between clusters.
#   Workflow: 24a prepare input (sc-core) -> in Docker: build Oracle -> PCA and KNN imputation -> fit GRN per cluster -> filter links, network scores -> ranking and comparison plots
# [Pipeline position] Upstream 05 -> this step (end point)
# [Environment] The 24a block runs in sc-core; the rest runs in Jupyter inside the Docker container kenjikamimoto126/celloracle_ubuntu:0.18.0
#   (see SKILL.md for how to start it; the container mounts the whole project folder, fair_utils.py included) (course ps30)
# [Input]
#   - results/05_annotate/objects/05_annotated.h5ad (layers['counts'], obs['cell_type'], obsm['X_umap'])
#   - Prior base GRN: bundled with CellOracle (human or mouse)
# [Output] All in results/24_celloracle/
#   - objects/: 24_input.h5ad (from 24a), 24_oracle.celloracle.oracle, 24_links.celloracle.links
#   - tables/: 24_grn_<cluster>.csv (filtered TF -> target links), 24_network_scores.csv
#   - figures/: 24_pca_cumulative_variance, 24_score_comparison_<cluster1>_vs_<cluster2>, 24_network_score_<TF>, plus ranking plots saved by CellOracle's own functions
# [How to read the results]
#   - degree_centrality_all: many connections = hub gene
#   - eigenvector_centrality: genes connected to other important nodes
#   - betweenness_centrality: genes acting as "bridges" for information flow
#   - score_comparison plot: TFs far from the diagonal differ strongly in importance between the two clusters, candidate state-specific regulators
#   - network_score_<TF>: how one TF's network importance changes across clusters
# [Notes]
#   - The course only covered network construction; for perturbation simulation (simulate_shift) see the official CellOracle tutorial
#   - Fitting the GRN for each cluster takes about half an hour

# %%
# 24a (sc-core environment)
# Step: prepare the CellOracle input: <= 30,000 cells, 3000 HVGs, raw counts in layers['raw_count']
# Input: results/05_annotate/objects/05_annotated.h5ad
# Output: results/24_celloracle/objects/24_input.h5ad
# import scanpy as sc
# import fair_utils as fair
# dirs = fair.setup_step("24_celloracle")
# adata = sc.read_h5ad("results/05_annotate/objects/05_annotated.h5ad")
# if adata.n_obs > 30000:
#     sc.pp.subsample(adata, n_obs=30000, random_state=0)
# sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3", layer="counts")
# adata = adata[:, adata.var.highly_variable].copy()
# adata.layers["raw_count"] = adata.layers["counts"].copy()
# adata.write(fair.out(dirs, "objects", "input.h5ad"))

# %%
# The following runs in Jupyter inside the Docker container (working directory = project folder)
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/24_celloracle/)
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import celloracle as co
import fair_utils as fair

IN_FILE = "results/24_celloracle/objects/24_input.h5ad"
CLUSTER_KEY = "cell_type"
EMBEDDING = "X_umap"
SPECIES = "mouse"                 # mouse: mouse scATAC atlas base GRN; human: human promoter base GRN
ALPHA = 10                        # GRN regularisation strength, larger = sparser
P_CUTOFF = 0.001
MAX_EDGES = 2000                  # maximum number of links kept per cluster
GOI = "Gata2"                     # transcription factor of interest
COMPARE = ("Stem", "Enterocyte")  # two clusters to compare centrality
dirs = fair.setup_step("24_celloracle")

# %%
# Step: read the data and the prior base GRN
# Input: IN_FILE
# Output: adata, base_grn
adata = sc.read_h5ad(IN_FILE)
if SPECIES == "mouse":
    base_grn = co.data.load_mouse_scATAC_atlas_base_GRN()
else:
    base_grn = co.data.load_human_promoter_base_GRN()
base_grn.head()

# %%
# Step: build the Oracle object, load raw counts and the prior network
# Input: adata, base_grn
# Output: oracle
adata.X = adata.layers["raw_count"].copy()
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name=CLUSTER_KEY, embedding_name=EMBEDDING)
oracle.import_TF_data(TF_info_matrix=base_grn)

# %%
# Step: PCA, choose the number of components at the cumulative variance elbow
# Input: oracle
# Output: n_comps; figures/24_pca_cumulative_variance.pdf/.png
oracle.perform_PCA()
cumulative = np.cumsum(oracle.pca.explained_variance_ratio_)[:100]
n_comps = min(np.where(np.diff(np.diff(cumulative) > 0.002))[0][0], 50)
plt.plot(cumulative)
plt.axvline(n_comps, color="black")
fair.save_fig(dirs, "pca_cumulative_variance")
n_comps

# %%
# Step: KNN imputation (k = 2.5% of the number of cells), save the oracle
# Input: oracle, n_comps
# Output: objects/24_oracle.celloracle.oracle
k = int(0.025 * oracle.adata.n_obs)
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k * 8, b_maxl=k * 4, n_jobs=4)
oracle.to_hdf5(str(fair.out(dirs, "objects", "oracle.celloracle.oracle")))

# %%
# Step: fit a GRN for each cluster (slow, about half an hour)
# Input: oracle
# Output: links (one link table per cluster)
links = oracle.get_links(cluster_name_for_GRN_unit=CLUSTER_KEY, alpha=ALPHA, verbose_level=10)

# %%
# Step: filter links, compute network scores, save links
# Input: links
# Output: objects/24_links.celloracle.links
links.filter_links(p=P_CUTOFF, weight="coef_abs", threshold_number=MAX_EDGES)
links.get_network_score()
links.to_hdf5(file_path=str(fair.out(dirs, "objects", "links.celloracle.links")))

# %%
# Step: save each cluster's link table and all network scores
# Input: links
# Output: tables/24_grn_<cluster>.csv, 24_network_scores.csv
for cluster, edges in links.filtered_links.items():
    edges.to_csv(fair.out(dirs, "tables", f"grn_{cluster}.csv"), index=False)
links.merged_score.to_csv(fair.out(dirs, "tables", "network_scores.csv"))

# %%
# Step: ranking plots of the most central genes per cluster (CellOracle's own function, saved straight to figures/)
# Input: links
# Output: ranking plots in figures/
for cluster in links.links_dict.keys():
    links.plot_scores_as_rank(cluster=cluster, n_gene=30, save=str(dirs["figures"]))

# %%
# Step: compare eigenvector centrality between two clusters
# Input: links, COMPARE
# Output: figures/24_score_comparison_<cluster1>_vs_<cluster2>.pdf/.png
links.plot_score_comparison_2D(value="eigenvector_centrality", cluster1=COMPARE[0], cluster2=COMPARE[1], percentile=98)
fair.save_fig(dirs, f"score_comparison_{COMPARE[0]}_vs_{COMPARE[1]}")

# %%
# Step: network score of one TF across clusters
# Input: links, GOI
# Output: figures/24_network_score_<TF>.pdf/.png
links.plot_score_per_cluster(goi=GOI)
fair.save_fig(dirs, f"network_score_{GOI}")

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
