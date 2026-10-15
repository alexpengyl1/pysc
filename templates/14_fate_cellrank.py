# 14 CellRank: terminal states, fate probabilities and lineage driver genes
#
# [Purpose]
#   Builds a cell-cell transition matrix from Palantir pseudotime (11) or RNA velocity (13),
#   finds stable macrostates and terminal states automatically, computes each cell's probability of reaching each terminal state, and finds the genes driving each lineage.
#   Workflow: transition matrix (directional kernel + similarity kernel) -> random walks to check direction -> GPCCA macrostates -> terminal states -> fate probabilities -> drivers -> gene trends
# [Pipeline position] Upstream 11 (Palantir) or 13 (scVelo) -> this step (end point)
# [Environment] Palantir pseudotime -> sc-palantir; RNA velocity -> sc-velo (course ps19)
# [Input]
#   - results/11_palantir/objects/11_palantir.h5ad (obs['palantir_pseudotime'])
#   - or results/13_scvelo/objects/13_velocity.h5ad (velocity layer; switch to the commented VelocityKernel block in the code)
# [Output] All in results/14_cellrank/
#   - objects/14_cellrank.h5ad
#   - tables/: 14_terminal_states.csv, 14_fate_probabilities.csv, 14_drivers_<lineage>.csv
#   - figures/: 14_random_walks, 14_macrostates_all, 14_macrostates_terminal, 14_fate_probabilities, 14_circular_projection, 14_gene_trends (CellRank's own plotting functions save PDF only)
# [How to read the results]
#   - random walks: simulated paths from the start cluster should reach the expected terminal clusters; if not, check the input pseudotime or velocity first
#   - macrostates: coarse-grained stable states; too small N_STATES merges lineages, too large splits off noise states
#   - fate probabilities: each cell's probability of reaching each terminal state
#   - circular projection: which fate each cell type leans towards
#   - lineage drivers: genes most correlated with a lineage's fate probability (high corr, small qval), candidate regulators

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/14_cellrank/)
import numpy as np
import pandas as pd
import scanpy as sc
import cellrank as cr
import fair_utils as fair

IN_FILE = "results/11_palantir/objects/11_palantir.h5ad"
TIME_KEY = "palantir_pseudotime"
CLUSTER_KEY = "cell_type"
START_CLUSTER = "Stem"
KERNEL_WEIGHT = 0.8        # weight of the directional kernel (pseudotime or velocity); the rest goes to the similarity kernel
N_STATES = 6
GENES = ["LGR5", "MKI67", "FABP2", "MUC2"]
SEED = 0
dirs = fair.setup_step("14_cellrank")
cr.settings.figdir = str(dirs["figures"])

# %%
# Step: read the data, recompute the neighbour graph
# Input: IN_FILE
# Output: adata
adata = sc.read_h5ad(IN_FILE)
sc.pp.neighbors(adata, random_state=SEED)
adata

# %%
# Step: build the transition matrix (directional kernel + expression similarity kernel)
# Input: adata.obs[TIME_KEY], adata.obsp['connectivities']
# Output: pk, kernel
pk = cr.kernels.PseudotimeKernel(adata, time_key=TIME_KEY).compute_transition_matrix()
ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
kernel = KERNEL_WEIGHT * pk + (1 - KERNEL_WEIGHT) * ck
# With RNA velocity use instead:
# vk = cr.kernels.VelocityKernel(adata).compute_transition_matrix()
# kernel = KERNEL_WEIGHT * vk + (1 - KERNEL_WEIGHT) * ck

# %%
# Step: simulate random walks to check the direction makes sense
# Input: pk
# Output: figures/14_random_walks.pdf
pk.plot_random_walks(n_sims=100, start_ixs={CLUSTER_KEY: START_CLUSTER}, basis="umap", seed=SEED,
                     save="14_random_walks.pdf")

# %%
# Step: GPCCA coarse-graining, plot all macrostates
# Input: kernel
# Output: g (estimator); figures/14_macrostates_all.pdf
g = cr.estimators.GPCCA(kernel)
g.fit(n_states=N_STATES, cluster_key=CLUSTER_KEY)
g.plot_macrostates(which="all", save="14_macrostates_all.pdf")

# %%
# Step: predict terminal states, save each cell's terminal state label
# Input: g
# Output: figures/14_macrostates_terminal.pdf; tables/14_terminal_states.csv
g.predict_terminal_states()
g.plot_macrostates(which="terminal", save="14_macrostates_terminal.pdf")
g.terminal_states.to_csv(fair.out(dirs, "tables", "terminal_states.csv"))

# %%
# Step: fate probabilities (plot + table)
# Input: g
# Output: figures/14_fate_probabilities.pdf; tables/14_fate_probabilities.csv
g.compute_fate_probabilities()
g.plot_fate_probabilities(same_plot=False, save="14_fate_probabilities.pdf")
fate = pd.DataFrame(np.asarray(g.fate_probabilities), index=adata.obs_names, columns=g.fate_probabilities.names)
fate.to_csv(fair.out(dirs, "tables", "fate_probabilities.csv"))

# %%
# Step: fate bias of each cell type (circular projection)
# Input: adata (fate probabilities)
# Output: figures/14_circular_projection.pdf
cr.pl.circular_projection(adata, keys=CLUSTER_KEY, legend_loc="right", save="14_circular_projection.pdf")

# %%
# Step: driver genes of each lineage
# Input: g
# Output: one tables/14_drivers_<lineage>.csv per lineage
for lineage in g.terminal_states.cat.categories:
    drivers = g.compute_lineage_drivers(lineages=[lineage])
    drivers.to_csv(fair.out(dirs, "tables", f"drivers_{lineage}.csv"))

# %%
# Step: gene trends (pyGAM, no R needed; switch to cr.models.GAMR for R mgcv)
# Input: adata, fate probabilities
# Output: figures/14_gene_trends.pdf
model = cr.models.GAM(adata)
cr.pl.gene_trends(adata, model=model, genes=GENES, time_key=TIME_KEY, same_plot=True, hide_cells=True,
                  save="14_gene_trends.pdf")

# %%
# Step: save
# Input: adata
# Output: objects/14_cellrank.h5ad
adata.write(fair.out(dirs, "objects", "cellrank.h5ad"))

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
