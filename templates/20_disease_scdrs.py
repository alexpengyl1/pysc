# 20 scDRS: find cells related to a disease / trait
#
# [Purpose]
#   Maps a disease (or trait) gene set from GWAS onto single-cell data, gives each cell a disease score,
#   compares it with many expression-matched random gene sets for significance, and tests which cell types are associated and whether there is heterogeneity within a cell type.
#   Workflow: normalise -> covariates -> read gene set -> regress covariates and bin genes -> score each trait -> UMAP -> cell-type-level tests
# [Pipeline position] Upstream 05 -> this step (end point)
# [Environment] sc-core (course ps26)
# [Input]
#   - results/05_annotate/objects/05_annotated.h5ad (layers['counts'])
#   - Gene set file .gs: two columns TRAIT, GENESET, genes formatted as "GENE1:weight,GENE2:weight";
#     made from GWAS summary statistics with MAGMA followed by scdrs munge-gs, or downloaded from the 74 traits provided by scDRS
# [Output] All in results/20_scdrs/, file names include the trait
#   - tables/: 20_score_<trait>.csv (each cell's score, p-value and 1000 control scores), 20_group_stats_<trait>.csv
#   - figures/20_umap_scores
# [How to read the results]
#   - norm_score: each cell's normalised disease score, higher = more associated; see where it concentrates on the UMAP
#   - group_stats table: small assoc_mcp = cell type significantly associated with the disease; small hetero_mcp = heterogeneous association within the type (possibly a subpopulation)
#   - Use an unrelated trait such as Height as a negative control to make sure the result is not universal
# [Notes]
#   - For mouse data use the ortholog conversion option of load_gs

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/20_scdrs/)
import pandas as pd
import scanpy as sc
import scdrs
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
GS_FILE = "data/geneset.gs"
SPECIES = "human"          # species of the expression data: human / mouse
GS_SPECIES = "human"       # species of the gene sets (MAGMA results are usually human)
N_CTRL = 1000
dirs = fair.setup_step("20_scdrs")

# %%
# Step: read the data, set X to raw counts, apply the filtering and normalisation scDRS recommends
# Input: IN_FILE
# Output: adata (X = log-normalised counts)
adata = sc.read_h5ad(IN_FILE)
adata.X = adata.layers["counts"].copy()
sc.pp.filter_cells(adata, min_genes=250)
sc.pp.filter_genes(adata, min_cells=50)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# %%
# Step: covariate table (intercept + number of genes detected per cell)
# Input: adata.obs
# Output: cov (cells x covariates)
cov = pd.DataFrame({"const": 1, "n_genes": adata.obs["n_genes"]}, index=adata.obs_names)

# %%
# Step: read the gene sets, convert orthologs if needed, intersect with expressed genes
# Input: GS_FILE
# Output: gene_sets (trait -> (gene list, weight list))
gene_sets = scdrs.util.load_gs(GS_FILE, src_species=GS_SPECIES, dst_species=SPECIES,
                               to_intersect=adata.var_names)
list(gene_sets.keys())

# %%
# Step: regress out covariates, bin genes by expression mean/variance
# Input: adata, cov
# Output: adata (modified in place)
scdrs.preprocess(adata, cov=cov, n_mean_bin=20, n_var_bin=20, copy=False)

# %%
# Step: score every cell for each trait
# Input: adata, gene_sets
# Output: scores (trait -> score table); one tables/20_score_<trait>.csv per trait
scores = {}
for trait, (genes, weights) in gene_sets.items():
    scores[trait] = scdrs.score_cell(adata, genes, gene_weight=weights, ctrl_match_key="mean_var",
                                     n_ctrl=N_CTRL, weight_opt="vs", return_ctrl_norm_score=True)
    scores[trait].to_csv(fair.out(dirs, "tables", f"score_{trait}.csv"))

# %%
# Step: scores on the UMAP
# Input: scores
# Output: figures/20_umap_scores.pdf/.png
for trait in scores:
    adata.obs[trait] = scores[trait]["norm_score"]
sc.pl.umap(adata, color=list(scores), cmap="RdBu_r", vmin=-5, vmax=5, show=False)
fair.save_fig(dirs, "umap_scores")

# %%
# Step: cell-type-level association and heterogeneity tests
# Input: adata, scores
# Output: one tables/20_group_stats_<trait>.csv per trait
for trait in scores:
    stats = scdrs.method.downstream_group_analysis(adata, df_full_score=scores[trait], group_cols=["cell_type"])
    stats["cell_type"].to_csv(fair.out(dirs, "tables", f"group_stats_{trait}.csv"))

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
