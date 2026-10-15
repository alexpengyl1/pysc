# 23 Metacells: aggregate similar cells into metacells
#
# [Purpose]
#   Uses the metacells algorithm from the Tanay lab to merge dozens of cells with almost identical transcriptional states into one "metacell".
#   Merged data are less noisy and smaller, which suits analyses sensitive to sparsity such as gene correlation and GRNs, and can be explored interactively in MCView.
#   Workflow: UMI distribution -> exclude noisy genes and low-quality cells -> mark lateral genes (cell cycle etc.) -> compute metacells -> transfer cell type labels -> UMAP -> save
# [Pipeline position] Upstream 05 -> this step (end point, or input for other analyses)
# [Environment] sc-metacell (course ps29)
# [Input]
#   - results/05_annotate/objects/05_annotated.h5ad (layers['counts']; all genes)
# [Output] All in results/23_metacells/
#   - objects/: 23_cells.h5ad (metacell of each cell in obs['metacell_name']), 23_metacells.h5ad (one row per metacell)
#   - tables/: 23_cell_to_metacell.csv, 23_metacell_summary.csv, 23_lateral_gene_modules.csv
#   - figures/: 23_total_umis_per_cell, 23_metacells_umap_cell_type
# [How to read the results]
#   - total_umis_per_cell: cells between the two threshold lines are kept
#   - Each metacell usually contains a few dozen cells; outlier cells are not assigned to any metacell
#   - Lateral genes (cell cycle, ribosomal, stress) do not drive grouping but are still used for outlier detection; check the lateral_gene_modules table and add whole modules to lateral if needed
#   - Metacell UMAP coloured by cell type: metacells of the same cell type should sit together
# [Notes]
#   - The course demo subsampled 50,000 cells; use all cells for real analyses, and lower the parallelism if memory runs out

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/23_metacells/)
import matplotlib.pyplot as plt
import scanpy as sc
import metacells as mc
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
MIN_UMIS = 800                 # course: 800-20000 UMIs
MAX_UMIS = 20000
MAX_EXCLUDED_FRACTION = 0.25   # maximum UMI fraction from excluded genes (mitochondrial etc.)
EXCLUDED_GENE_NAMES = ["XIST", "MALAT1"]
EXCLUDED_GENE_PATTERNS = ["MT-.*"]              # mouse "mt-.*"
LATERAL_GENE_NAMES = ["AURKA", "MCM3", "MCM4", "MCM7", "MKI67", "PCNA", "RRM2", "SMC4", "TPX2",
                      "FOS", "HSP90AB1", "TXN", "DUSP1"]
LATERAL_GENE_PATTERNS = ["RP[LS].*"]
SEED = 123456
dirs = fair.setup_step("23_metacells")

# %%
# Step: read the data, keep only raw counts (metacells needs float32 UMIs)
# Input: IN_FILE
# Output: full
annotated = sc.read_h5ad(IN_FILE)
full = sc.AnnData(annotated.layers["counts"].astype("float32"), obs=annotated.obs, var=annotated.var[[]])
mc.ut.set_name(full, "full")
mc.ut.top_level(full)
full

# %%
# Step: distribution of total UMIs per cell with threshold lines
# Input: full
# Output: figures/23_total_umis_per_cell.pdf/.png
plt.hist(full.X.sum(axis=1).A1, bins=100, log=True)
plt.axvline(MIN_UMIS, color="green")
plt.axvline(MAX_UMIS, color="red")
plt.xlabel("UMIs per cell")
fair.save_fig(dirs, "total_umis_per_cell")

# %%
# Step: exclude specified genes and let metacells exclude noisy genes automatically
# Input: full
# Output: full.var['excluded_gene']
mc.pl.exclude_genes(full, excluded_gene_names=EXCLUDED_GENE_NAMES,
                    excluded_gene_patterns=EXCLUDED_GENE_PATTERNS, random_seed=SEED)

# %%
# Step: mark cells to exclude by total UMIs and excluded-gene fraction
# Input: full
# Output: full.obs['excluded_cell']
mc.tl.compute_excluded_gene_umis(full)
mc.pl.exclude_cells(full, properly_sampled_min_cell_total=MIN_UMIS,
                    properly_sampled_max_cell_total=MAX_UMIS,
                    properly_sampled_max_excluded_genes_fraction=MAX_EXCLUDED_FRACTION)

# %%
# Step: extract the clean data
# Input: full
# Output: cells
cells = mc.pl.extract_clean_data(full, name="cells")
mc.ut.top_level(cells)
cells

# %%
# Step: mark lateral genes
# Input: cells
# Output: cells.var['lateral_gene']
mc.pl.mark_lateral_genes(cells, lateral_gene_names=LATERAL_GENE_NAMES,
                         lateral_gene_patterns=LATERAL_GENE_PATTERNS)

# %%
# Step: find gene modules correlated with lateral genes, save the table for manual review (whole modules can be added to lateral, see course ps29)
# Input: cells
# Output: cells.var['lateral_genes_module']; tables/23_lateral_gene_modules.csv
mc.pl.relate_to_lateral_genes(cells, random_seed=SEED)
cells.var[["lateral_gene", "lateral_genes_module"]].to_csv(fair.out(dirs, "tables", "lateral_gene_modules.csv"))

# %%
# Step: compute metacells with the divide-and-conquer pipeline (lower the parallelism by hand if memory runs out)
# Input: cells
# Output: cells.obs['metacell_name']
mc.pl.set_max_parallel_piles(mc.pl.guess_max_parallel_piles(cells))
with mc.ut.progress_bar():
    mc.pl.divide_and_conquer_pipeline(cells, random_seed=SEED)

# %%
# Step: collect metacells and transfer cell type labels to them (majority)
# Input: cells
# Output: metacells (metacells x genes)
metacells = mc.pl.collect_metacells(cells, name="metacells", random_seed=SEED)
mc.tl.convey_obs_to_group(adata=cells, gdata=metacells, property_name="cell_type", method=mc.ut.most_frequent)
metacells

# %%
# Step: save the cell -> metacell table and the metacell summary
# Input: cells, metacells
# Output: tables/23_cell_to_metacell.csv, 23_metacell_summary.csv
cells.obs[["cell_type", "metacell_name"]].to_csv(fair.out(dirs, "tables", "cell_to_metacell.csv"))
metacells.obs.to_csv(fair.out(dirs, "tables", "metacell_summary.csv"))

# %%
# Step: compute metacell UMAP and other display data (also read by MCView), colour by cell type
# Input: cells, metacells
# Output: metacells.obs['x'], ['y']; figures/23_metacells_umap_cell_type.pdf/.png
mc.pl.compute_for_mcview(adata=cells, gdata=metacells, random_seed=SEED)
sc.pl.scatter(metacells, x="x", y="y", color="cell_type", show=False)
fair.save_fig(dirs, "metacells_umap_cell_type")

# %%
# Step: save
# Input: cells, metacells
# Output: objects/23_cells.h5ad, 23_metacells.h5ad
cells.write_h5ad(fair.out(dirs, "objects", "cells.h5ad"))
metacells.write_h5ad(fair.out(dirs, "objects", "metacells.h5ad"))

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
