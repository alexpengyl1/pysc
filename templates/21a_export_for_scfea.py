# 21a Export the expression matrix for scFEA
#
# [Purpose]
#   scFEA is a separate command-line tool (its own scFEA environment) that needs a genes x cells csv expression matrix.
#   This step exports the log-normalised expression in that format from sc-core, together with each cell's cell type; cells can be subsampled when there are many.
#   (Course ps27 exported the data slot from Seurat; this is the scanpy equivalent)
# [Pipeline position] Upstream 05 -> this step -> 21_metabolism_scfea.sh
# [Environment] sc-core
# [Input]
#   - results/05_annotate/objects/05_annotated.h5ad (layers['log1p'], obs['cell_type'])
# [Output] All in results/21_scfea/ (shared by 21a and 21)
#   - objects/21_expr.csv: genes x cells (input to 21_metabolism_scfea.sh)
#   - tables/21_cell_info.csv: cell -> cell_type (for summarising fluxes by cell type later)
#   - logs/run_info_21a.json, pip_freeze_21a.txt
# [How to read the results]
#   - This step only converts formats; check the printed matrix size for the expected numbers of genes and cells
# [Notes]
#   - scFEA trains a model on every run, so more cells = slower; adjust N_CELLS to your computer

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/21_scfea/)
import scanpy as sc
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
N_CELLS = 2000     # scFEA trains a model on every run, so subsample first
SEED = 0
dirs = fair.setup_step("21_scfea")

# %%
# Step: read the data and subsample
# Input: IN_FILE
# Output: adata (N_CELLS cells)
adata = sc.read_h5ad(IN_FILE)
sc.pp.subsample(adata, n_obs=min(N_CELLS, adata.n_obs), random_state=SEED)

# %%
# Step: export the genes x cells log-normalised matrix and cell annotation
# Input: adata
# Output: objects/21_expr.csv, tables/21_cell_info.csv
adata.to_df(layer="log1p").T.to_csv(fair.out(dirs, "objects", "expr.csv"))
adata.obs[["cell_type"]].to_csv(fair.out(dirs, "tables", "cell_info.csv"))

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info_21a.json, logs/pip_freeze_21a.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()}, script="21a")
