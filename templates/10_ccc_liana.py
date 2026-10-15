# 10 Cell-cell communication: LIANA+ multi-method consensus ranking
#
# [Purpose]
#   Runs several cell-cell communication methods at once (CellPhoneDB, CellChat, NATMI, Connectome, log2FC etc.)
#   and combines them into one consensus ranking to reduce the bias of any single method, then filters specific, strong ligand-receptor interactions and plots them.
# [Pipeline position] Upstream 05 -> this step (end point; can be cross-checked with 09 CellPhoneDB)
# [Environment] liana (course ps15)
# [Input]
#   - results/05_annotate/objects/05_annotated.h5ad (X = log1p; obs['cell_type'])
#   - SOURCES / TARGETS at the top: the sender and receiver cells to focus on
# [Output] All in results/10_liana/
#   - tables/: 10_liana_rank_aggregate.csv (all results), 10_liana_top_specific.csv (specificity_rank <= 0.05)
#   - figures/: 10_dotplot_sources_targets, 10_circle_plot_counts
# [How to read the results]
#   - magnitude_rank: consensus rank of interaction strength, smaller = stronger
#   - specificity_rank: how specific the interaction is to this pair of cell types, smaller = more specific
#   - Usually filter by specificity_rank <= 0.05 first, then sort by magnitude_rank
#   - Network plot: thicker line = more specific interactions between two cell types
# [Notes]
#   - Use RESOURCE = "consensus" for human and "mouseconsensus" for mouse

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/10_liana/)
import scanpy as sc
import liana as li
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
CELLTYPE_KEY = "cell_type"
RESOURCE = "consensus"            # mouse: "mouseconsensus"
EXPR_PROP = 0.1
SPECIFICITY_CUTOFF = 0.05
N_PERMS = 1000
SEED = 0
SOURCES = ["Fibroblast", "Myeloid"]
TARGETS = ["T", "Epithelial"]
dirs = fair.setup_step("10_liana")

# %%
# Step: read the data
# Input: IN_FILE
# Output: adata
adata = sc.read_h5ad(IN_FILE)
adata

# %%
# Step: multi-method consensus ranking (X is log1p, so use_raw=False)
# Input: adata
# Output: adata.uns['liana_res']; tables/10_liana_rank_aggregate.csv
li.mt.rank_aggregate(adata, groupby=CELLTYPE_KEY, resource_name=RESOURCE, expr_prop=EXPR_PROP,
                     use_raw=False, n_perms=N_PERMS, seed=SEED, verbose=True)
liana_res = adata.uns["liana_res"]
liana_res.to_csv(fair.out(dirs, "tables", "liana_rank_aggregate.csv"), index=False)
liana_res.head()

# %%
# Step: keep specific interactions only, sorted by strength
# Input: liana_res
# Output: tables/10_liana_top_specific.csv
top = liana_res[liana_res["specificity_rank"] <= SPECIFICITY_CUTOFF].sort_values("magnitude_rank")
top.to_csv(fair.out(dirs, "tables", "liana_top_specific.csv"), index=False)
top.head()

# %%
# Step: dot plot of interactions for chosen source / target cells
# Input: adata.uns['liana_res']
# Output: figures/10_dotplot_sources_targets.pdf/.png
dotplot = li.pl.dotplot(adata=adata, colour="magnitude_rank", size="specificity_rank",
                        inverse_colour=True, inverse_size=True,
                        source_labels=SOURCES, target_labels=TARGETS,
                        filter_fun=lambda x: x["specificity_rank"] <= SPECIFICITY_CUTOFF,
                        top_n=20, orderby="magnitude_rank", orderby_ascending=True, figure_size=(10, 7))
fair.save_fig(dirs, "dotplot_sources_targets", fig=dotplot)

# %%
# Step: network plot of specific interaction counts between cell types
# Input: adata.uns['liana_res']
# Output: figures/10_circle_plot_counts.pdf/.png
li.pl.circle_plot(adata, groupby=CELLTYPE_KEY, score_key="magnitude_rank", inverse_score=True,
                  source_labels=SOURCES, filter_fun=lambda x: x["specificity_rank"] <= SPECIFICITY_CUTOFF,
                  pivot_mode="counts", figure_size=(8, 8))
fair.save_fig(dirs, "circle_plot_counts")

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
