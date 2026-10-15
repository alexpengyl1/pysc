# 07 Within one cell type: differential expression between conditions + GSEA + Enrichr
#
# [Purpose]
#   In one chosen cell type, compares two conditions (CASE vs CONTROL) to find differentially expressed genes,
#   then looks for changed pathways in two ways: GSEA prerank (using the ranking of all genes) and Enrichr (using only significant up / down genes).
#   Workflow: download gene set library -> subset cell type -> Wilcoxon DE -> GSEA prerank -> Enrichr (up and down separately) -> plots
# [Pipeline position] Upstream 05 -> this step (end point; change CELLTYPE to repeat for each cell type)
# [Environment] gseapy (course ps12); needs internet to download gene sets
# [Input]
#   - results/05_annotate/objects/05_annotated.h5ad (X = log1p; obs['cell_type'], obs['condition'])
#   - CELLTYPE, CASE, CONTROL and the gene set library name (GO / KEGG / MSigDB) at the top
# [Output] All in results/07_gseapy/, file names include the cell type (e.g. 07_deg_cd14_mono.csv)
#   - tables/: 07_deg_<celltype>.csv (DE results for all genes), 07_prerank_<celltype>.csv, 07_enrichr_up_<celltype>.csv, 07_enrichr_down_<celltype>.csv
#   - figures/: 07_deg_top25_<celltype>, 07_prerank_top5_<celltype>, 07_enrichr_up_dotplot_<celltype>, 07_enrichr_down_dotplot_<celltype>
# [How to read the results]
#   - DE table: scores > 0 = up in CASE; look at pvals_adj and logfoldchanges
#   - prerank: NES > 0 = pathway up overall in CASE, NES < 0 = down; FDR q-val < 0.25 is the usual GSEA cut-off
#   - Enrichment curve: a peak on the left = pathway genes concentrated at the up-regulated end
#   - Enrichr dot plots: pathways with Adjusted P-value < 0.05; read up and down separately
# [Notes]
#   - Per-cell tests treat cells as replicates, so p-values come out too small; with several samples, pseudobulk per sample is safer

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/07_gseapy/)
import scanpy as sc
import gseapy as gp
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
CELLTYPE = "CD14 Mono"
CONDITION_KEY = "condition"
CASE, CONTROL = "STIM", "CTRL"
ORGANISM = "Human"                         # or "Mouse"
LIBRARY = "GO_Biological_Process_2023"      # gp.get_library_name(organism=ORGANISM) lists them all
SEED = 0
TAG = CELLTYPE.lower().replace(" ", "_")   # cell type label used in file names
dirs = fair.setup_step("07_gseapy")

# %%
# Step: download the gene set library from Enrichr
# Input: LIBRARY, ORGANISM
# Output: gene_sets (dict: pathway name -> gene list)
gene_sets = gp.get_library(name=LIBRARY, organism=ORGANISM)
len(gene_sets)

# %%
# Step: read the data, take one cell type
# Input: IN_FILE
# Output: bdata (subset of that cell type); prints the number of cells per group
adata = sc.read_h5ad(IN_FILE)
bdata = adata[adata.obs["cell_type"] == CELLTYPE].copy()
bdata.obs[CONDITION_KEY].value_counts()

# %%
# Step: CASE vs CONTROL differential expression (Wilcoxon)
# Input: bdata
# Output: deg (DE table); tables/07_deg_<celltype>.csv
sc.tl.rank_genes_groups(bdata, groupby=CONDITION_KEY, groups=[CASE], reference=CONTROL,
                        method="wilcoxon", use_raw=False)
deg = sc.get.rank_genes_groups_df(bdata, group=CASE)
deg.to_csv(fair.out(dirs, "tables", f"deg_{TAG}.csv"), index=False)
deg.head()

# %%
# Step: plot of the top 25 DE genes
# Input: bdata.uns['rank_genes_groups']
# Output: figures/07_deg_top25_<celltype>.pdf/.png
sc.pl.rank_genes_groups(bdata, n_genes=25, sharey=False, show=False)
fair.save_fig(dirs, f"deg_top25_{TAG}")

# %%
# Step: GSEA prerank on genes ranked by Wilcoxon score
# Input: deg (names, scores)
# Output: prerank (result object); tables/07_prerank_<celltype>.csv
prerank = gp.prerank(rnk=deg[["names", "scores"]], gene_sets=gene_sets,
                     permutation_num=1000, seed=SEED, outdir=None)
prerank.res2d.to_csv(fair.out(dirs, "tables", f"prerank_{TAG}.csv"), index=False)
prerank.res2d.head(10)

# %%
# Step: enrichment curves of the top 5 pathways
# Input: prerank
# Output: figures/07_prerank_top5_<celltype>.pdf/.png
prerank.plot(terms=prerank.res2d.Term[:5])
fair.save_fig(dirs, f"prerank_top5_{TAG}")

# %%
# Step: Enrichr on significant up / down genes separately (background = all detected genes)
# Input: deg
# Output: enr_up, enr_down; tables/07_enrichr_up_<celltype>.csv, 07_enrichr_down_<celltype>.csv
sig = deg[deg.pvals_adj < 0.05]
background = bdata.var_names.tolist()
enr_up = gp.enrichr(gene_list=sig[sig.logfoldchanges > 0].names.tolist(), gene_sets=gene_sets,
                    background=background, outdir=None)
enr_down = gp.enrichr(gene_list=sig[sig.logfoldchanges < 0].names.tolist(), gene_sets=gene_sets,
                      background=background, outdir=None)
enr_up.res2d.to_csv(fair.out(dirs, "tables", f"enrichr_up_{TAG}.csv"), index=False)
enr_down.res2d.to_csv(fair.out(dirs, "tables", f"enrichr_down_{TAG}.csv"), index=False)

# %%
# Step: Enrichr dot plot for up-regulated genes
# Input: enr_up
# Output: figures/07_enrichr_up_dotplot_<celltype>.pdf/.png
gp.dotplot(enr_up.res2d, title=f"{CELLTYPE} up in {CASE}", figsize=(4, 6))
fair.save_fig(dirs, f"enrichr_up_dotplot_{TAG}")

# %%
# Step: Enrichr dot plot for down-regulated genes
# Input: enr_down
# Output: figures/07_enrichr_down_dotplot_<celltype>.pdf/.png
gp.dotplot(enr_down.res2d, title=f"{CELLTYPE} down in {CASE}", figsize=(4, 6))
fair.save_fig(dirs, f"enrichr_down_dotplot_{TAG}")

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
