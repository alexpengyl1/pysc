# 13 RNA velocity (scVelo)
#
# [Purpose]
#   Uses the ratio of unspliced to spliced mRNA to infer which direction each cell's gene expression is moving in (RNA velocity),
#   predicting the direction of cell state transitions, and estimates an internal "latent time" and the genes driving these changes.
#   Workflow: read and merge loom files -> splicing proportions -> normalise and neighbour moments -> dynamical model -> velocity stream plot -> confidence -> latent time -> driver genes
# [Pipeline position] Upstream 05 (+ loom files from velocyto) -> this step -> downstream 14_fate_cellrank (optional)
# [Environment] sc-velo (course ps18)
# [Input]
#   - One loom per sample (spliced / unspliced counts), generated beforehand with velocyto from the Cell Ranger output (needs BAM, genome GTF, optional repeat-mask GTF):
#     pip install velocyto (needs samtools), then velocyto run10x -m repeat_msk.gtf <cellranger_sample_dir> genes.gtf;
#     the file appears in <cellranger_sample_dir>/velocyto/<sample>.loom
#   - results/05_annotate/objects/05_annotated.h5ad
# [Output] All in results/13_scvelo/
#   - objects/13_velocity.h5ad (usable by 14)
#   - tables/: 13_spliced_proportions.csv, 13_velocity_confidence.csv, 13_top_likelihood_genes.csv, 13_dynamical_genes_cell_type.csv
#   - figures/: 13_proportions_cell_type, 13_velocity_stream_cell_type, 13_velocity_confidence_length, 13_latent_time, 13_heatmap_top300_latent_time, 13_phase_portraits_top10
# [How to read the results]
#   - proportions: unspliced is usually about 10-25%; much lower makes velocity unreliable
#   - Stream plot: arrows = inferred direction of state change; check against single-gene phase portraits, do not rely on arrows alone
#   - velocity_confidence: directions in low-confidence regions are not trustworthy
#   - latent_time: an "internal time" based on transcriptional dynamics
#   - Driver genes: genes with high fit_likelihood have the clearest dynamics; a clean arc in the phase portrait (spliced vs unspliced) = good fit
# [Notes]
#   - Terminally differentiated or quiescent populations have weak velocity, and their arrows are easily misleading

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/13_scvelo/)
import anndata as ad
import scanpy as sc
import scvelo as scv
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
LOOMS = {"S1": "data/velocyto/S1.loom", "S2": "data/velocyto/S2.loom"}
MIN_SHARED_COUNTS = 20
N_TOP_GENES = 2000
MODE = "dynamical"      # if short of time use "stochastic" and delete the recover_dynamics and latent_time blocks
N_JOBS = 4
scv.settings.verbosity = 3
dirs = fair.setup_step("13_scvelo")

# %%
# Step: read the annotated data and loom files; rename loom cells to match 05 (<barcode>-1-<sample>)
# Input: IN_FILE, LOOMS
# Output: adata, ldata (loom files of all samples merged)
adata = sc.read_h5ad(IN_FILE)
looms = {}
for name, path in LOOMS.items():
    looms[name] = scv.read(path, cache=True)
    looms[name].var_names_make_unique()
    looms[name].obs_names = [bc.split(":")[1].replace("x", "-1") + "-" + name for bc in looms[name].obs_names]
ldata = ad.concat(looms)
ldata

# %%
# Step: merge spliced/unspliced into the annotated data
# Input: adata, ldata
# Output: adata (adds layers['spliced'], ['unspliced'])
adata = scv.utils.merge(adata, ldata)
adata

# %%
# Step: spliced / unspliced proportions per cell type (plot + table)
# Input: adata
# Output: figures/13_proportions_cell_type.pdf/.png; tables/13_spliced_proportions.csv
scv.pl.proportions(adata, groupby="cell_type", show=False)
fair.save_fig(dirs, "proportions_cell_type")
proportions = adata.obs[["cell_type"]].assign(
    spliced=adata.layers["spliced"].sum(axis=1).A1, unspliced=adata.layers["unspliced"].sum(axis=1).A1)
proportions.groupby("cell_type")[["spliced", "unspliced"]].sum().to_csv(fair.out(dirs, "tables", "spliced_proportions.csv"))

# %%
# Step: filter and normalise, compute neighbour moments
# Input: adata
# Output: adata.layers['Ms'], ['Mu']
scv.pp.filter_and_normalize(adata, min_shared_counts=MIN_SHARED_COUNTS, n_top_genes=N_TOP_GENES)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# %%
# Step: estimate velocity with the dynamical model (slow)
# Input: adata
# Output: adata.layers['velocity'], uns['velocity_graph']
scv.tl.recover_dynamics(adata, n_jobs=N_JOBS)
scv.tl.velocity(adata, mode=MODE)
scv.tl.velocity_graph(adata)

# %%
# Step: velocity stream plot
# Input: adata (obsm['X_umap'])
# Output: figures/13_velocity_stream_cell_type.pdf/.png
scv.pl.velocity_embedding_stream(adata, basis="umap", color="cell_type", show=False)
fair.save_fig(dirs, "velocity_stream_cell_type")

# %%
# Step: velocity confidence and length (plot + table)
# Input: adata
# Output: figures/13_velocity_confidence_length.pdf/.png; tables/13_velocity_confidence.csv
scv.tl.velocity_confidence(adata)
scv.pl.scatter(adata, c=["velocity_length", "velocity_confidence"], cmap="coolwarm", perc=[5, 95], show=False)
fair.save_fig(dirs, "velocity_confidence_length")
adata.obs[["cell_type", "velocity_length", "velocity_confidence"]].to_csv(fair.out(dirs, "tables", "velocity_confidence.csv"))

# %%
# Step: latent time
# Input: adata
# Output: obs['latent_time']; figures/13_latent_time.pdf/.png
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color="latent_time", color_map="gnuplot", show=False)
fair.save_fig(dirs, "latent_time")

# %%
# Step: driver genes ranked by fit likelihood (table)
# Input: adata.var['fit_likelihood']
# Output: top_genes; tables/13_top_likelihood_genes.csv
top_genes = adata.var["fit_likelihood"].sort_values(ascending=False)
top_genes.to_csv(fair.out(dirs, "tables", "top_likelihood_genes.csv"))
top_genes = top_genes.index

# %%
# Step: heatmap of the top 300 driver genes along latent time
# Input: adata, top_genes
# Output: figures/13_heatmap_top300_latent_time.pdf/.png
scv.pl.heatmap(adata, var_names=top_genes[:300], sortby="latent_time", col_color="cell_type", n_convolve=100, show=False)
fair.save_fig(dirs, "heatmap_top300_latent_time")

# %%
# Step: phase portraits of the top 10 driver genes
# Input: adata, top_genes
# Output: figures/13_phase_portraits_top10.pdf/.png
scv.pl.scatter(adata, basis=top_genes[:10], ncols=5, color="cell_type", frameon=False, show=False)
fair.save_fig(dirs, "phase_portraits_top10")

# %%
# Step: cell-type-specific dynamical genes
# Input: adata
# Output: tables/13_dynamical_genes_cell_type.csv
scv.tl.rank_dynamical_genes(adata, groupby="cell_type")
scv.get_df(adata, "rank_dynamical_genes/names").to_csv(fair.out(dirs, "tables", "dynamical_genes_cell_type.csv"))

# %%
# Step: save
# Input: adata
# Output: objects/13_velocity.h5ad
adata.write(fair.out(dirs, "objects", "velocity.h5ad"))

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
