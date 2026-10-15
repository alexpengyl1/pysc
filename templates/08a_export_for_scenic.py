# 08a Export the expression matrix for pySCENIC
#
# [Purpose]
#   pySCENIC only runs in a Python 3.7 environment and cannot read new h5ad files, so this step exports the raw counts to csv from sc-core,
#   together with each cell's cell type, for use in 08b and 08c. Genes with very low expression are dropped, and cells can be subsampled to a cap.
# [Pipeline position] Upstream 05 -> this step -> 08b_pyscenic -> 08c_scenic_plot
# [Environment] sc-core (course ps13)
# [Input]
#   - results/05_annotate/objects/05_annotated.h5ad (layers['counts'], obs['cell_type'])
# [Output] All in results/08_scenic/ (shared by 08a/08b/08c)
#   - objects/08_expr_counts.csv: cells x genes, raw counts (input to 08b)
#   - tables/08_cell_info.csv: cell -> cell_type (used in 08c)
#   - logs/run_info_08a.json, pip_freeze_08a.txt
# [How to read the results]
#   - This step only converts formats; check the printed matrix size for the expected numbers of cells and genes
# [Notes]
#   - The MAX_CELLS cap is only for speed; use all cells for real analyses when possible (the course demo used only 100 cells, which is not enough for conclusions)

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/08_scenic/)
import scanpy as sc
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
MIN_CELL_FRACTION = 0.01   # gene must be expressed in at least 1% of cells
MAX_CELLS = 20000          # GRNBoost2 run time grows with cell number; subsample above this
SEED = 0
dirs = fair.setup_step("08_scenic")

# %%
# Step: read the data, drop very lowly expressed genes, subsample if needed
# Input: IN_FILE
# Output: adata
adata = sc.read_h5ad(IN_FILE)
sc.pp.filter_genes(adata, min_cells=int(MIN_CELL_FRACTION * adata.n_obs))
if adata.n_obs > MAX_CELLS:
    sc.pp.subsample(adata, n_obs=MAX_CELLS, random_state=SEED)
adata

# %%
# Step: export the raw count matrix and cell annotation
# Input: adata.layers['counts'], adata.obs
# Output: objects/08_expr_counts.csv, tables/08_cell_info.csv
adata.to_df(layer="counts").to_csv(fair.out(dirs, "objects", "expr_counts.csv"))
adata.obs[["cell_type"]].to_csv(fair.out(dirs, "tables", "cell_info.csv"))

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info_08a.json, logs/pip_freeze_08a.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()}, script="08a")
