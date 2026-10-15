# 09 Cell-cell communication: CellPhoneDB v5 (statistical method) + ktplotspy plots
#
# [Purpose]
#   Uses ligand and receptor expression across cell types to find likely signalling between cell types (ligand-receptor pairs),
#   tests which interactions are specific with permutations, and counts significant interactions between each pair of cell types.
#   Workflow: download database -> (mouse) convert to human orthologs -> write input files -> permutation test -> count significant interactions -> heatmap and dot plot
# [Pipeline position] Upstream 05 -> this step (end point; can be cross-checked with 10 LIANA)
# [Environment] cpdb (course ps14)
# [Input]
#   - results/05_annotate/objects/05_annotated.h5ad (X = log1p, not scaled; obs['cell_type'])
#   - Optional: data/microenvironment.tsv (cell_type, microenvironment) to compare only cells in the same microenvironment; data/active_TFs.tsv (cell_type, TF)
# [Output] All in results/09_cellphonedb/
#   - objects/: 09_counts_log1p.h5ad, 09_metadata.tsv (CellPhoneDB input files)
#   - tables/: raw CellPhoneDB output statistical_analysis_pvalues / means / significant_means / interaction_scores / deconvoluted (file names carry 09), 09_n_significant_by_pair.csv
#   - figures/: 09_heatmap_n_significant, 09_dotplot_<source>_to_<target>
# [How to read the results]
#   - Heatmap: number of significant interactions between each pair of cell types; more = more communication (only "likely", not experimental evidence)
#   - significant_means: rows = ligand|receptor, columns = cellA|cellB; a value = significant, empty = not significant
#   - interaction_scores: 0-100, higher = more specific
#   - Dot plot: dot size = significance, colour = mean expression; pick pathways relevant to your question for validation
# [Notes]
#   - The database has human genes only. Convert mouse data to human gene names with an ortholog table (Ensembl BioMart export: mouse_symbol, human_symbol); do not just uppercase the names

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/09_cellphonedb/)
import pandas as pd
import anndata as ad
import ktplotspy as kpy
from cellphonedb.utils import db_utils
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
DB_DIR = "data/cpdb_db"
DB_VERSION = "v5.0.0"
CELLTYPE_KEY = "cell_type"
SPECIES = "human"              # for mouse use "mouse" and prepare the ortholog table ORTHOLOG_FILE
ORTHOLOG_FILE = "data/mouse_human_orthologs.csv"
ITERATIONS = 1000
EXPR_THRESHOLD = 0.1           # genes expressed in < 10% of a cell type count as not expressed
PVALUE = 0.05
SEED = 42
SOURCE = "Fibroblast"                       # for plotting: ligand source cells
TARGETS = "Myeloid|T"                       # for plotting: receptor cells, separated by |
GENES = ["TGFB1", "CXCL12", "CSF1"]
dirs = fair.setup_step("09_cellphonedb")

# %%
# Step: download the CellPhoneDB database (once only)
# Input: none
# Output: data/cpdb_db/cellphonedb.zip
db_utils.download_database(DB_DIR, DB_VERSION)

# %%
# Step: read the annotated data
# Input: IN_FILE
# Output: adata
adata = ad.read_h5ad(IN_FILE)
adata

# %%
# Step (mouse only): convert gene names to human orthologs, drop genes without an ortholog or with duplicates
# Input: adata, ORTHOLOG_FILE
# Output: adata (var_names = human gene names)
if SPECIES == "mouse":
    orthologs = pd.read_csv(ORTHOLOG_FILE).drop_duplicates("mouse_symbol").set_index("mouse_symbol")
    adata = adata[:, adata.var_names.isin(orthologs.index)].copy()
    adata.var_names = orthologs.loc[adata.var_names, "human_symbol"].values
    adata = adata[:, ~adata.var_names.duplicated()].copy()

# %%
# Step: write the expression matrix and cell annotation CellPhoneDB needs
# Input: adata
# Output: objects/09_counts_log1p.h5ad, objects/09_metadata.tsv
counts_file = fair.out(dirs, "objects", "counts_log1p.h5ad")
meta_file = fair.out(dirs, "objects", "metadata.tsv")
adata.write(counts_file)
pd.DataFrame({"barcode_sample": adata.obs_names, "cell_type": adata.obs[CELLTYPE_KEY].values}).to_csv(
    meta_file, sep="\t", index=False)

# %%
# Step: run the statistical analysis (permutation test)
# Input: database zip, metadata, counts_log1p
# Output: cpdb_results (dict); CellPhoneDB result tables written to tables/
cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path=f"{DB_DIR}/cellphonedb.zip",
    meta_file_path=str(meta_file),
    counts_file_path=str(counts_file),
    counts_data="hgnc_symbol",
    # microenvs_file_path="data/microenvironment.tsv",   # optional: pair only within the same microenvironment
    # active_tfs_file_path="data/active_TFs.tsv",        # optional: CellSign
    score_interactions=True,
    iterations=ITERATIONS,
    threshold=EXPR_THRESHOLD,
    threads=4,
    debug_seed=SEED,
    pvalue=PVALUE,
    output_path=str(dirs["tables"]),
    output_suffix="09",
)
list(cpdb_results.keys())

# %%
# Step: number of significant interactions for each pair of cell types
# Input: cpdb_results['pvalues']
# Output: tables/09_n_significant_by_pair.csv
pvals = cpdb_results["pvalues"].set_index("interacting_pair").filter(like="|")
n_sig = (pvals < PVALUE).sum().rename("n_significant").sort_values(ascending=False)
n_sig.to_csv(fair.out(dirs, "tables", "n_significant_by_pair.csv"))
n_sig.head()

# %%
# Step: heatmap of significant interaction counts between cell types
# Input: cpdb_results['pvalues']
# Output: figures/09_heatmap_n_significant.pdf/.png
grid = kpy.plot_cpdb_heatmap(pvals=cpdb_results["pvalues"], figsize=(6, 6), title="Number of significant interactions")
fair.save_fig(dirs, "heatmap_n_significant", fig=grid)

# %%
# Step: dot plot of interactions for chosen cell types and genes
# Input: adata, cpdb_results
# Output: figures/09_dotplot_<source>_to_<target>.pdf/.png
dotplot = kpy.plot_cpdb(
    adata=adata,
    cell_type1=SOURCE,
    cell_type2=TARGETS,
    means=cpdb_results["means"],
    pvals=cpdb_results["pvalues"],
    celltype_key=CELLTYPE_KEY,
    genes=GENES,
    interaction_scores=cpdb_results["interaction_scores"],
    figsize=(10, 5),
)
fair.save_fig(dirs, f"dotplot_{SOURCE}_to_{TARGETS.replace('|', '_')}".lower(), fig=dotplot)

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
