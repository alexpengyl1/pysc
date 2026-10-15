# 16 MAGIC expression imputation
#
# [Purpose]
#   Many genes read 0 in single-cell data (dropout), which hides relationships between genes.
#   MAGIC smooths selected genes using information from similar cells so that gene-gene relationships (e.g. VIM up, CDH1 down in EMT) become visible.
# [Pipeline position] Upstream 05 -> this step (end point)
# [Environment] sc-core (course ps21)
# [Input]
#   - results/05_annotate/objects/05_annotated.h5ad (X = log1p)
#   - GENES at the top: genes to impute and show
# [Output] All in results/16_magic/
#   - tables/16_magic_imputed_<genes>.csv: imputed expression of GENES only (cells x genes)
#   - figures/: 16_scatter_before_magic, 16_scatter_after_magic
# [How to read the results]
#   - Before imputation: many points sit on the axes (zeros) and no relationship is visible
#   - After imputation: a continuous trend should appear; if there is still no relationship, the two genes are simply not related
# [Notes]
#   - Imputed values are only for looking at gene relationships, trends and figures, never for differential expression or statistical tests

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/16_magic/)
import scanpy as sc
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
GENES = ["VIM", "CDH1", "ZEB1"]    # first two on x / y, third for colour
KNN = 5        # smaller = more local, less smoothing
T = 3          # diffusion steps, larger = smoother
SEED = 0
TAG = "_".join(GENES).lower()
dirs = fair.setup_step("16_magic")

# %%
# Step: read the data
# Input: IN_FILE
# Output: adata
adata = sc.read_h5ad(IN_FILE)

# %%
# Step: MAGIC-impute the chosen genes and save the imputed values
# Input: adata.X
# Output: adata_magic (cells x GENES); tables/16_magic_imputed_<genes>.csv
adata_magic = sc.external.pp.magic(adata, name_list=GENES, knn=KNN, t=T, random_state=SEED, copy=True)
adata_magic.to_df().to_csv(fair.out(dirs, "tables", f"magic_imputed_{TAG}.csv"))

# %%
# Step: gene relationship before imputation
# Input: adata
# Output: figures/16_scatter_before_magic.pdf/.png
sc.pl.scatter(adata, x=GENES[0], y=GENES[1], color=GENES[2], title="Before MAGIC", show=False)
fair.save_fig(dirs, "scatter_before_magic")

# %%
# Step: gene relationship after imputation
# Input: adata_magic
# Output: figures/16_scatter_after_magic.pdf/.png
sc.pl.scatter(adata_magic, x=GENES[0], y=GENES[1], color=GENES[2], title="After MAGIC", show=False)
fair.save_fig(dirs, "scatter_after_magic")

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
