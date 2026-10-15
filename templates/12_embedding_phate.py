# 12 PHATE embedding
#
# [Purpose]
#   Draws cells in two dimensions with PHATE. Compared with UMAP, PHATE preserves continuous differentiation or time-course structure (branches, transitions) better.
#   It is a visualisation only and produces no pseudotime or statistics.
# [Pipeline position] Upstream 05 -> this step (end point)
# [Environment] sc-core (course ps17)
# [Input]
#   - results/05_annotate/objects/05_annotated.h5ad (X = log1p)
# [Output] All in results/12_phate/
#   - tables/12_phate_coordinates.csv: PHATE1 / PHATE2 coordinates of each cell
#   - figures/12_phate_cell_type_sample
# [How to read the results]
#   - Colour by cell type / time point / sample and look for continuous trajectories and branches; branch tips are usually terminal cell types
#   - Parameters: smaller KNN -> finer local structure; for large datasets (> 100,000 cells) KNN can go up to 30-40; T controls smoothing, larger = smoother
# [Notes]
#   - Distances and shapes in the plot are for visual display only; support conclusions with pseudotime (11) or RNA velocity (13)

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/12_phate/)
import pandas as pd
import scanpy as sc
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
KNN = 5
T = "auto"
COLOR_BY = ["cell_type", "sample"]
SEED = 0
dirs = fair.setup_step("12_phate")

# %%
# Step: read the data
# Input: IN_FILE
# Output: adata
adata = sc.read_h5ad(IN_FILE)

# %%
# Step: compute the PHATE embedding and save coordinates
# Input: adata.X
# Output: adata.obsm['X_phate']; tables/12_phate_coordinates.csv
sc.external.tl.phate(adata, k=KNN, t=T, n_components=2, n_jobs=-2, random_state=SEED)
pd.DataFrame(adata.obsm["X_phate"], index=adata.obs_names, columns=["PHATE1", "PHATE2"]).to_csv(
    fair.out(dirs, "tables", "phate_coordinates.csv"))

# %%
# Step: plot PHATE
# Input: adata.obsm['X_phate']
# Output: figures/12_phate_cell_type_sample.pdf/.png
sc.external.pl.phate(adata, color=COLOR_BY, show=False)
fair.save_fig(dirs, "phate_cell_type_sample")

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
