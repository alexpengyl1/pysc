# 11 Palantir: pseudotime, differentiation potential and terminal fate probabilities
#
# [Purpose]
#   For one differentiation lineage (e.g. stem cells -> several mature cell types), infers each cell's progress along differentiation (pseudotime),
#   how open its fate still is (entropy), and its probability of reaching each terminal state, and plots gene trends along each branch.
#   Workflow: subset the lineage -> diffusion map -> MAGIC imputation -> check start and terminal cells -> run Palantir -> branch trajectories -> gene trends
# [Pipeline position] Upstream 05 -> this step -> downstream 14_fate_cellrank (optional)
# [Environment] sc-palantir (course ps16)
# [Input]
#   - The lineage subset of results/05_annotate/objects/05_annotated.h5ad (obs['cell_type'], obsm['X_umap'])
#   - START_CELL at the top: barcode of the start cell (e.g. the cell with the highest stemness markers in the stem/progenitor cluster)
#   - TERMINAL_STATES at the top: barcodes of representative cells for each terminal state
# [Output] All in results/11_palantir/
#   - objects/11_palantir.h5ad: obs has palantir_pseudotime, palantir_entropy; obsm has palantir_fate_probabilities; layers['MAGIC_imputed_data'] (used by 14)
#   - tables/: 11_pseudotime_entropy.csv, 11_fate_probabilities.csv
#   - figures/: 11_start_terminal_cells, 11_palantir_results, 11_branch_selection, 11_trajectories, 11_gene_trends, 11_gene_trend_heatmaps
# [How to read the results]
#   - start_terminal_cells plot: first make sure the start and terminal cells sit in the right clusters; wrong choices make everything after wrong
#   - pseudotime: progress along differentiation from the start cell, 0 -> 1
#   - entropy: high = multipotent, fate not yet decided; should fall along differentiation
#   - fate probabilities: each cell's probability of reaching each terminal state, summing to 1
#   - gene trends: how genes change with pseudotime along each branch; look for branch-specific up / down genes
# [Notes]
#   - Include only cells of one lineage; unrelated cell types distort the trajectory
#   - MAGIC-imputed values are only for drawing trends, not for statistical tests

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/11_palantir/)
import pandas as pd
import scanpy as sc
import palantir
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
LINEAGE_TYPES = ["Stem", "TA", "Enterocyte", "Goblet"]   # cell types of the lineage to analyse
START_CELL = "AAACCCAAGXXXXXXX-1-S1"                     # barcode of the start cell
TERMINAL_STATES = {                                       # representative cell barcode: terminal state name
    "AAACGAAXXXXXXX-1-S1": "Enterocyte",
    "AAAGTCCXXXXXXX-1-S2": "Goblet",
}
N_WAYPOINTS = 500
GENES = ["LGR5", "MKI67", "FABP2", "MUC2"]
SEED = 0
dirs = fair.setup_step("11_palantir")

# %%
# Step: read the data, take the lineage subset, recompute PCA
# Input: IN_FILE
# Output: adata (lineage subset)
adata = sc.read_h5ad(IN_FILE)
adata = adata[adata.obs["cell_type"].isin(LINEAGE_TYPES)].copy()
sc.pp.pca(adata, mask_var="highly_variable", random_state=SEED)
adata

# %%
# Step: diffusion map and multiscale space
# Input: adata.obsm['X_pca']
# Output: adata.obsm['DM_EigenVectors_multiscaled'] etc.
palantir.utils.run_diffusion_maps(adata, n_components=5)
palantir.utils.determine_multiscale_space(adata)

# %%
# Step: MAGIC imputation (only for gene trends and display)
# Input: adata
# Output: adata.layers['MAGIC_imputed_data']
palantir.utils.run_magic_imputation(adata)

# %%
# Step: check the positions of the start and terminal cells on the UMAP
# Input: START_CELL, TERMINAL_STATES
# Output: figures/11_start_terminal_cells.pdf/.png
terminal_states = pd.Series(TERMINAL_STATES)
palantir.plot.highlight_cells_on_umap(adata, pd.concat([terminal_states, pd.Series({START_CELL: "start"})]))
fair.save_fig(dirs, "start_terminal_cells")

# %%
# Step: run Palantir
# Input: adata, START_CELL, terminal_states
# Output: adata.obs['palantir_pseudotime'], ['palantir_entropy'], obsm['palantir_fate_probabilities']
palantir.core.run_palantir(adata, START_CELL, num_waypoints=N_WAYPOINTS, terminal_states=terminal_states, seed=SEED)

# %%
# Step: overview of pseudotime, entropy and fate probabilities
# Input: adata
# Output: figures/11_palantir_results.pdf/.png
palantir.plot.plot_palantir_results(adata, s=3)
fair.save_fig(dirs, "palantir_results")

# %%
# Step: save pseudotime, entropy and fate probabilities per cell
# Input: adata
# Output: tables/11_pseudotime_entropy.csv, 11_fate_probabilities.csv
adata.obs[["cell_type", "palantir_pseudotime", "palantir_entropy"]].to_csv(fair.out(dirs, "tables", "pseudotime_entropy.csv"))
pd.DataFrame(adata.obsm["palantir_fate_probabilities"], index=adata.obs_names).to_csv(fair.out(dirs, "tables", "fate_probabilities.csv"))

# %%
# Step: select representative cells of each branch
# Input: adata
# Output: adata.obsm['branch_masks']; figures/11_branch_selection.pdf/.png
palantir.presults.select_branch_cells(adata, q=0.01, eps=0.01)
palantir.plot.plot_branch_selection(adata)
fair.save_fig(dirs, "branch_selection")

# %%
# Step: trajectory of each branch
# Input: adata
# Output: figures/11_trajectories.pdf/.png
palantir.plot.plot_trajectories(adata, pseudotime_interval=(0, 0.9))
fair.save_fig(dirs, "trajectories")

# %%
# Step: gene expression trends along each branch (line plots)
# Input: adata.layers['MAGIC_imputed_data']
# Output: figures/11_gene_trends.pdf/.png
palantir.presults.compute_gene_trends(adata, expression_key="MAGIC_imputed_data")
palantir.plot.plot_gene_trends(adata, GENES)
fair.save_fig(dirs, "gene_trends")

# %%
# Step: gene trend heatmaps
# Input: adata
# Output: figures/11_gene_trend_heatmaps.pdf/.png
palantir.plot.plot_gene_trend_heatmaps(adata, GENES)
fair.save_fig(dirs, "gene_trend_heatmaps")

# %%
# Step: save
# Input: adata
# Output: objects/11_palantir.h5ad
adata.write(fair.out(dirs, "objects", "palantir.h5ad"))

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
