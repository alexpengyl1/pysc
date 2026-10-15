# 22 muon: 10x Multiome (RNA + ATAC) analysis and integration
#
# [Purpose]
#   When gene expression (RNA) and chromatin accessibility (ATAC) are measured in the same cells, first QC, normalise, reduce and cluster each modality separately,
#   then integrate them into a joint representation with MOFA+ (shared and modality-specific factors) and WNN (per-cell weighted combination of both neighbour graphs).
#   Workflow: RNA QC/clustering/markers -> ATAC QC (fragment length, nucleosome signal, TSS enrichment)/clustering -> intersect cells -> MOFA+ -> WNN UMAP
# [Pipeline position] Standalone (reads Cell Ranger ARC output directly, does not depend on 01-05)
# [Environment] sc-multi (course ps28)
# [Input]
#   - filtered_feature_bc_matrix.h5 from Cell Ranger ARC, with atac_fragments.tsv.gz and its .tbi index in the same folder
# [Output] All in results/22_muon/
#   - objects/22_multiome.h5mu: rna and atac modalities; obsm['X_mofa'], WNN UMAP; objects/22_mofa_model.hdf5
#   - tables/: 22_qc_cells_kept.csv, 22_markers_rna_leiden.csv
#   - figures/: 22_rna_qc_violin, 22_rna_umap_leiden, 22_atac_qc_histogram, 22_atac_fragment_length, 22_atac_nucleosome_signal, 22_atac_tss_enrichment, 22_atac_umap_leiden, 22_wnn_umap
# [How to read the results]
#   - RNA QC: as in 02
#   - Fragment length plot: should show a nucleosome-free peak < 147 bp and a mono-nucleosome peak around 200 bp (periodicity)
#   - nucleosome_signal: cells with very high values are poor quality and are removed
#   - TSS enrichment: the sharper the peak at transcription start sites, the better the data
#   - MOFA factors: capture variation shared by both modalities as well as modality-specific variation
#   - WNN UMAP: colour by RNA clusters and by ATAC clusters; agreement between them means the integration is reasonable
# [Notes]
#   - Mouse mitochondrial genes start with mt-

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/22_muon/)
import pandas as pd
import scanpy as sc
import muon as mu
from muon import atac as ac
import fair_utils as fair

H5_FILE = "data/multiome/filtered_feature_bc_matrix.h5"
MT_PREFIX = "MT-"
RNA_MIN_GENES, RNA_MAX_GENES, RNA_MAX_PCT_MT = 200, 8000, 10
ATAC_MIN_COUNTS, ATAC_MAX_COUNTS = 1000, 80000
N_FACTORS = 20
RESOLUTION = 0.5
SEED = 0
dirs = fair.setup_step("22_muon")

# %%
# Step: read the multiome data
# Input: H5_FILE (the fragments file is linked automatically)
# Output: mdata (modalities rna, atac)
mdata = mu.read_10x_h5(H5_FILE)
mdata.var_names_make_unique()
cells_before = {"rna": mdata["rna"].n_obs, "atac": mdata["atac"].n_obs}
mdata

# %%
# Step: RNA QC metrics and violin plots
# Input: mdata['rna']
# Output: rna; figures/22_rna_qc_violin.pdf/.png
rna = mdata.mod["rna"]
rna.var["mt"] = rna.var_names.str.startswith(MT_PREFIX)
sc.pp.calculate_qc_metrics(rna, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(rna, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], multi_panel=True, show=False)
fair.save_fig(dirs, "rna_qc_violin")

# %%
# Step: filter RNA cells
# Input: rna
# Output: rna (filtered)
mu.pp.filter_obs(rna, "n_genes_by_counts", lambda x: (x >= RNA_MIN_GENES) & (x < RNA_MAX_GENES))
mu.pp.filter_obs(rna, "pct_counts_mt", lambda x: x < RNA_MAX_PCT_MT)

# %%
# Step: RNA normalisation, HVGs, PCA, clustering, UMAP
# Input: rna
# Output: rna.obs['leiden'], obsm['X_pca'], ['X_umap']; figures/22_rna_umap_leiden.pdf/.png
rna.layers["counts"] = rna.X.copy()
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
rna.layers["log1p"] = rna.X.copy()
sc.pp.highly_variable_genes(rna, min_mean=0.02, max_mean=4, min_disp=0.5)
sc.pp.scale(rna, max_value=10)
sc.tl.pca(rna, random_state=SEED)
sc.pp.neighbors(rna, n_pcs=20, random_state=SEED)
sc.tl.leiden(rna, resolution=RESOLUTION, flavor="igraph", n_iterations=2, directed=False, random_state=SEED)
sc.tl.umap(rna, random_state=SEED)
sc.pl.umap(rna, color="leiden", legend_loc="on data", show=False)
fair.save_fig(dirs, "rna_umap_leiden")

# %%
# Step: RNA markers per cluster (using the log1p layer)
# Input: rna
# Output: tables/22_markers_rna_leiden.csv
sc.tl.rank_genes_groups(rna, "leiden", method="wilcoxon", layer="log1p", use_raw=False)
sc.get.rank_genes_groups_df(rna, group=None).to_csv(fair.out(dirs, "tables", "markers_rna_leiden.csv"), index=False)

# %%
# Step: distributions of ATAC QC metrics
# Input: mdata['atac']
# Output: atac; figures/22_atac_qc_histogram.pdf/.png
atac = mdata.mod["atac"]
sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)
mu.pl.histogram(atac, ["n_genes_by_counts", "total_counts"])
fair.save_fig(dirs, "atac_qc_histogram")

# %%
# Step: filter ATAC peaks and cells
# Input: atac
# Output: atac (filtered)
mu.pp.filter_var(atac, "n_cells_by_counts", lambda x: x >= 10)
mu.pp.filter_obs(atac, "total_counts", lambda x: (x >= ATAC_MIN_COUNTS) & (x <= ATAC_MAX_COUNTS))

# %%
# Step: fragment length distribution
# Input: atac (fragments file)
# Output: figures/22_atac_fragment_length.pdf/.png
ac.pl.fragment_histogram(atac, region="chr1:1-2000000")
fair.save_fig(dirs, "atac_fragment_length")

# %%
# Step: nucleosome signal
# Input: atac
# Output: atac.obs['nucleosome_signal']; figures/22_atac_nucleosome_signal.pdf/.png
ac.tl.nucleosome_signal(atac, n=1e6)
mu.pl.histogram(atac, "nucleosome_signal")
fair.save_fig(dirs, "atac_nucleosome_signal")

# %%
# Step: TSS enrichment
# Input: mdata
# Output: atac.obs['tss_score']; figures/22_atac_tss_enrichment.pdf/.png
tss = ac.tl.tss_enrichment(mdata, n_tss=1000)
ac.pl.tss_enrichment(tss)
fair.save_fig(dirs, "atac_tss_enrichment")

# %%
# Step: ATAC normalisation, highly variable peaks, PCA, clustering, UMAP
# Input: atac
# Output: atac.obs['leiden'], obsm['X_pca'], ['X_umap']; figures/22_atac_umap_leiden.pdf/.png
atac.layers["counts"] = atac.X.copy()
sc.pp.normalize_total(atac, target_sum=1e4)
sc.pp.log1p(atac)
sc.pp.highly_variable_genes(atac, min_mean=0.05, max_mean=1.5, min_disp=0.5)
sc.pp.scale(atac, max_value=10)
sc.tl.pca(atac, random_state=SEED)
sc.pp.neighbors(atac, n_pcs=20, random_state=SEED)
sc.tl.leiden(atac, resolution=RESOLUTION, flavor="igraph", n_iterations=2, directed=False, random_state=SEED)
sc.tl.umap(atac, random_state=SEED)
sc.pl.umap(atac, color="leiden", legend_loc="on data", show=False)
fair.save_fig(dirs, "atac_umap_leiden")

# %%
# Step: intersect the two modalities, keep cells present in both, record cell numbers
# Input: mdata
# Output: mdata (cells aligned); tables/22_qc_cells_kept.csv
mdata.update()
mu.pp.intersect_obs(mdata)
pd.DataFrame({"before": cells_before, "after": {"rna": mdata["rna"].n_obs, "atac": mdata["atac"].n_obs}}).to_csv(
    fair.out(dirs, "tables", "qc_cells_kept.csv"))

# %%
# Step: MOFA+ factor integration
# Input: mdata
# Output: mdata.obsm['X_mofa']; objects/22_mofa_model.hdf5
mu.tl.mofa(mdata, n_factors=N_FACTORS, outfile=str(fair.out(dirs, "objects", "mofa_model.hdf5")))

# %%
# Step: WNN weighted neighbour graph and UMAP (neighbours recomputed for each modality first)
# Input: mdata
# Output: mdata.obsm['X_umap']; figures/22_wnn_umap.pdf/.png
sc.pp.neighbors(mdata["rna"])
sc.pp.neighbors(mdata["atac"])
mu.pp.neighbors(mdata, key_added="wnn")
mu.tl.umap(mdata, neighbors_key="wnn", random_state=SEED)
mu.pl.umap(mdata, color=["rna:leiden", "atac:leiden"], show=False)
fair.save_fig(dirs, "wnn_umap")

# %%
# Step: save
# Input: mdata
# Output: objects/22_multiome.h5mu
mdata.write(fair.out(dirs, "objects", "multiome.h5mu"))

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
