# 04 scVI / scANVI integration + scib-metrics benchmark
#
# [Purpose]
#   Uses the deep generative model scVI to learn a batch-corrected low-dimensional representation (X_scVI) from raw counts, then clusters and draws a UMAP.
#   When some cell labels already exist, scANVI can refine the model; with labels, scib-metrics can also compare integration results quantitatively.
#   Workflow: take HVGs -> train scVI -> check training curve -> cluster/UMAP on the latent space -> (optional) scANVI -> (optional) scib benchmark -> save
# [Pipeline position] Upstream 02 -> this step -> downstream 05_annotate (alternative to 03; point IN_FILE in 05 to this output)
# [Environment] sc-scvi (PyTorch; course ps10)
# [Input]
#   - results/02_qc_cluster/objects/02_processed.h5ad: layers['counts'] (must be raw counts), obs['sample']; scANVI / scib also need obs['cell_type']
# [Output] All in results/04_scvi/
#   - objects/04_scvi.h5ad: obsm['X_scVI'] (optionally X_scANVI), leiden, X_umap
#   - objects/04_scvi_model/: the trained model, can be reloaded without retraining
#   - figures/: 04_training_elbo, 04_umap_scvi_sample_leiden, (optional) 04_scib_results
#   - tables/: 04_cells_per_cluster_sample.csv, (optional) 04_scib_metrics.csv
# [How to read the results]
#   - Training curve: ELBO should fall smoothly and level off; a validation curve that rises again = overfitting, lower max_epochs in train()
#   - UMAP: mixed by sample and separated by cell type = good integration
#   - scib table: high Batch correction = batches removed well; high Bio conservation = biological differences kept; read them together, with X_pca as the unintegrated baseline
# [Notes]
#   - Runs without a GPU, but slowly when there are many cells

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/04_scvi/)
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import scvi
import fair_utils as fair

IN_FILE = "results/02_qc_cluster/objects/02_processed.h5ad"
BATCH_KEY = "sample"
LABEL_KEY = "cell_type"      # only needed for scANVI / scib
N_LATENT = 30
RESOLUTION = 0.8
SEED = 0
scvi.settings.seed = SEED
dirs = fair.setup_step("04_scvi")

# %%
# Step: read the data, keep only highly variable genes for modelling
# Input: IN_FILE
# Output: adata (all genes), adata_hvg (HVG subset)
adata = sc.read_h5ad(IN_FILE)
adata_hvg = adata[:, adata.var.highly_variable].copy()
adata_hvg

# %%
# Step: register the data and train scVI (input must be raw counts)
# Input: adata_hvg.layers['counts']
# Output: model (trained SCVI)
scvi.model.SCVI.setup_anndata(adata_hvg, layer="counts", batch_key=BATCH_KEY)
model = scvi.model.SCVI(adata_hvg, n_layers=2, n_latent=N_LATENT, gene_likelihood="nb")
model.train()

# %%
# Step: training curve
# Input: model.history
# Output: figures/04_training_elbo.pdf/.png
model.history["elbo_train"].plot()
fair.save_fig(dirs, "training_elbo")

# %%
# Step: get the latent space, cluster and UMAP
# Input: model
# Output: adata.obsm['X_scVI'], obs['leiden'], obsm['X_umap']; figures/04_umap_scvi_sample_leiden.pdf/.png
adata.obsm["X_scVI"] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI", random_state=SEED)
sc.tl.leiden(adata, resolution=RESOLUTION, flavor="igraph", n_iterations=2, directed=False, random_state=SEED)
sc.tl.umap(adata, min_dist=0.3, random_state=SEED)
sc.pl.umap(adata, color=[BATCH_KEY, "leiden"], show=False)
fair.save_fig(dirs, "umap_scvi_sample_leiden")

# %%
# Step: cells from each sample in each cluster
# Input: adata.obs
# Output: tables/04_cells_per_cluster_sample.csv
pd.crosstab(adata.obs["leiden"], adata.obs[BATCH_KEY]).to_csv(fair.out(dirs, "tables", "cells_per_cluster_sample.csv"))

# %%
# Optional scANVI: fine-tune on top of scVI when some or all cell type labels are known
# Step: initialise scANVI from scVI and train
# Input: model, adata_hvg.obs[LABEL_KEY] (unannotated cells set to "Unknown")
# Output: adata.obsm['X_scANVI']
# adata_hvg.obs[LABEL_KEY] = adata.obs[LABEL_KEY]
# scanvi_model = scvi.model.SCANVI.from_scvi_model(model, adata=adata_hvg, labels_key=LABEL_KEY, unlabeled_category="Unknown")
# scanvi_model.train(max_epochs=20, n_samples_per_label=100)
# adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata_hvg)

# %%
# Optional scib benchmark: needs cell type labels
# Step: compare integration quality of X_pca, X_scVI, X_scANVI
# Input: adata (obs[BATCH_KEY], obs[LABEL_KEY], each obsm)
# Output: tables/04_scib_metrics.csv; figures/04_scib_results.pdf/.png
# from scib_metrics.benchmark import Benchmarker
# bm = Benchmarker(adata, batch_key=BATCH_KEY, label_key=LABEL_KEY,
#                  embedding_obsm_keys=["X_pca", "X_scVI", "X_scANVI"], n_jobs=-1)
# bm.benchmark()
# bm.get_results(min_max_scale=False).to_csv(fair.out(dirs, "tables", "scib_metrics.csv"))
# bm.plot_results_table(min_max_scale=False)
# fair.save_fig(dirs, "scib_results")

# %%
# Step: save the model and the data
# Input: model, adata
# Output: objects/04_scvi_model/, objects/04_scvi.h5ad
model.save(fair.out(dirs, "objects", "scvi_model"), overwrite=True)
adata.write(fair.out(dirs, "objects", "scvi.h5ad"))

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
