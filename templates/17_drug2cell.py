# 17 drug2cell: find cell-type-specific drug target signals
#
# [Purpose]
#   Scores every cell for the target gene set of each drug in ChEMBL (mean expression of the targets)
#   to find drugs whose targets are specifically high in a cell type, pointing to the cells a drug may act on, or to possible off-target effects.
#   Workflow: drug target scores -> drugs specific to each cell type in score space -> dot plot -> hypergeometric test of cell type markers vs target sets
# [Pipeline position] Upstream 05 -> this step (end point)
# [Environment] sc-core (course ps22)
# [Input]
#   - results/05_annotate/objects/05_annotated.h5ad (X = log1p; obs['cell_type'])
#   - ChEMBL target database: bundled with drug2cell
# [Output] All in results/17_drug2cell/
#   - objects/17_drug_scores.h5ad: cells x drugs score matrix
#   - tables/: 17_drug_markers_cell_type.csv, 17_hypergeom_<celltype>.csv
#   - figures/17_drug_dotplot_cell_type
# [How to read the results]
#   - Dot plot: a high score for a drug in a cell type -> the drug may act on these cells (including possible off-target effects)
#   - drug_markers table: drugs specific to each cell type; prefer high scores and small pvals_adj
#   - hypergeom table: whether the cell type's marker genes overlap the drug's targets significantly, a second line of evidence
# [Notes]
#   - The target database uses human gene names; convert mouse data to human orthologs first (see the ortholog block in 09_ccc_cellphonedb.py)
#   - Results are expression-based hints and need experimental or literature support

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/17_drug2cell/)
import scanpy as sc
import drug2cell as d2c
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
CELLTYPE_KEY = "cell_type"
N_TOP = 5
dirs = fair.setup_step("17_drug2cell")

# %%
# Step: read the data
# Input: IN_FILE
# Output: adata
adata = sc.read_h5ad(IN_FILE)

# %%
# Step: score every cell with the ChEMBL drug target sets, save the score matrix
# Input: adata.X (log1p)
# Output: drug_scores (cells x drugs AnnData); objects/17_drug_scores.h5ad
d2c.score(adata, use_raw=False)
drug_scores = adata.uns.pop("drug2cell")
drug_scores.write(fair.out(dirs, "objects", "drug_scores.h5ad"))
drug_scores

# %%
# Step: drugs specific to each cell type in drug score space
# Input: drug_scores
# Output: tables/17_drug_markers_cell_type.csv
sc.tl.rank_genes_groups(drug_scores, groupby=CELLTYPE_KEY, method="wilcoxon")
sc.get.rank_genes_groups_df(drug_scores, group=None).to_csv(fair.out(dirs, "tables", "drug_markers_cell_type.csv"), index=False)

# %%
# Step: dot plot of the top N_TOP specific drugs per cell type
# Input: drug_scores
# Output: figures/17_drug_dotplot_cell_type.pdf/.png
sc.pl.rank_genes_groups_dotplot(drug_scores, n_genes=N_TOP, swap_axes=True, dendrogram=False, show=False)
fair.save_fig(dirs, "drug_dotplot_cell_type")

# %%
# Step: markers in gene expression space, for over-representation analysis
# Input: adata
# Output: adata.uns['rank_genes_groups']
sc.tl.rank_genes_groups(adata, groupby=CELLTYPE_KEY, method="wilcoxon", use_raw=False)

# %%
# Step: hypergeometric test of markers vs drug target sets
# Input: adata
# Output: one tables/17_hypergeom_<celltype>.csv per cell type
overrep = d2c.hypergeometric(adata)
for celltype, table in overrep.items():
    table.to_csv(fair.out(dirs, "tables", f"hypergeom_{celltype.replace(' ', '_')}.csv"))

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
