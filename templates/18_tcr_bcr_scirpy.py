# 18 Scirpy: TCR / BCR immune receptor repertoire analysis
#
# [Purpose]
#   Combines the T/B cell receptor sequences from 10x VDJ with the transcriptome, defines clonotypes,
#   and analyses which T/B cells are clonally expanded, clonal diversity per cell type, V(D)J gene usage, clones shared between samples, and networks of similar clones.
#   Workflow: read VDJ -> build MuData with the transcriptome -> chain pairing QC -> drop multichain -> clonotypes and expansion -> diversity, V(D)J, CDR3, shared clones -> TCRdist network
# [Pipeline position] Upstream 05 (+ Cell Ranger vdj output) -> this step (end point)
# [Environment] sc-multi (course ps23)
# [Input]
#   - results/05_annotate/objects/05_annotated.h5ad (containing T or B cells)
#   - Cell Ranger vdj filtered_contig_annotations.csv for each sample (paths in VDJ_FILES, sample names as in 01)
# [Output] All in results/18_scirpy/
#   - objects/18_scirpy.h5mu
#   - tables/: 18_chain_pairing_sample.csv, 18_clonotypes_per_cell.csv, 18_clone_sizes.csv, 18_repertoire_overlap.csv
#   - figures/: 18_chain_pairing_sample, 18_clonal_expansion_cell_type, 18_alpha_diversity_cell_type, 18_umap_clonal_expansion, 18_vdj_usage, 18_spectratype_cell_type, 18_repertoire_overlap_sample, 18_clonotype_network_tcrdist
# [How to read the results]
#   - chain_pairing: multichain cells are mostly doublets and are removed; orphan = only one chain detected; single pair is the ideal case
#   - clonal_expansion: clone size classes; clusters rich in expanded clones are usually effector or tissue-resident cells
#   - alpha_diversity (normalised Shannon entropy): low = dominated by a few clones (clonal expansion)
#   - repertoire_overlap: clones shared between samples or clusters, suggesting a common antigen response
#   - TCRdist network: clones with similar amino acid sequences cluster together and may recognise the same antigen

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/18_scirpy/)
import pandas as pd
import anndata as ad
import scanpy as sc
import muon as mu
import scirpy as ir
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
VDJ_FILES = {
    "S1": "data/S1/filtered_contig_annotations.csv",
    "S2": "data/S2/filtered_contig_annotations.csv",
}
TCRDIST_CUTOFF = 15
dirs = fair.setup_step("18_scirpy")

# %%
# Step: read the transcriptome and each sample's VDJ, cell names in the same format as 05 (<barcode>-<sample>)
# Input: IN_FILE, VDJ_FILES
# Output: adata, airr (merged receptor data)
adata = sc.read_h5ad(IN_FILE)
airrs = {}
for name, path in VDJ_FILES.items():
    airrs[name] = ir.io.read_10x_vdj(path)
airr = ad.concat(airrs, index_unique="-")
airr

# %%
# Step: build the MuData
# Input: adata, airr
# Output: mdata (modalities gex, airr)
mdata = mu.MuData({"gex": adata, "airr": airr})
mdata

# %%
# Step: chain indexing and chain pairing QC (plot + table)
# Input: mdata
# Output: figures/18_chain_pairing_sample.pdf/.png; tables/18_chain_pairing_sample.csv
ir.pp.index_chains(mdata)
ir.tl.chain_qc(mdata)
ir.pl.group_abundance(mdata, groupby="airr:chain_pairing", target_col="gex:sample")
fair.save_fig(dirs, "chain_pairing_sample")
pd.crosstab(mdata.obs["airr:chain_pairing"], mdata.obs["gex:sample"]).to_csv(fair.out(dirs, "tables", "chain_pairing_sample.csv"))

# %%
# Step: remove multichain cells
# Input: mdata
# Output: mdata (filtered)
mu.pp.filter_obs(mdata, "airr:chain_pairing", lambda x: x != "multichain")

# %%
# Step: define clonotypes by identical CDR3 nucleotide sequences, compute clonal expansion
# Input: mdata
# Output: obs['airr:clone_id'], ['airr:clonal_expansion']
ir.pp.ir_dist(mdata)
ir.tl.define_clonotypes(mdata, receptor_arms="all", dual_ir="primary_only")
ir.tl.clonal_expansion(mdata)

# %%
# Step: save per-cell clone information and the size of each clone
# Input: mdata.obs
# Output: tables/18_clonotypes_per_cell.csv, 18_clone_sizes.csv
clone_cols = ["gex:sample", "gex:cell_type", "airr:clone_id", "airr:clone_id_size", "airr:clonal_expansion"]
mdata.obs[clone_cols].to_csv(fair.out(dirs, "tables", "clonotypes_per_cell.csv"))
mdata.obs.groupby(["airr:clone_id", "gex:sample"]).size().rename("n_cells").to_csv(fair.out(dirs, "tables", "clone_sizes.csv"))

# %%
# Step: clonal expansion by cell type
# Input: mdata
# Output: figures/18_clonal_expansion_cell_type.pdf/.png
ir.pl.clonal_expansion(mdata, target_col="clone_id", groupby="gex:cell_type", breakpoints=(1, 2, 5))
fair.save_fig(dirs, "clonal_expansion_cell_type")

# %%
# Step: clonal diversity by cell type
# Input: mdata
# Output: figures/18_alpha_diversity_cell_type.pdf/.png
ir.pl.alpha_diversity(mdata, metric="normalized_shannon_entropy", groupby="gex:cell_type")
fair.save_fig(dirs, "alpha_diversity_cell_type")

# %%
# Step: clonal expansion on the UMAP
# Input: mdata
# Output: figures/18_umap_clonal_expansion.pdf/.png
mu.pl.embedding(mdata, basis="gex:umap", color="airr:clonal_expansion", show=False)
fair.save_fig(dirs, "umap_clonal_expansion")

# %%
# Step: V(D)J gene usage
# Input: mdata
# Output: figures/18_vdj_usage.pdf/.png
ir.pl.vdj_usage(mdata, full_combination=False, max_ribbons=30)
fair.save_fig(dirs, "vdj_usage")

# %%
# Step: CDR3 length distribution
# Input: mdata
# Output: figures/18_spectratype_cell_type.pdf/.png
ir.pl.spectratype(mdata, color="gex:cell_type", viztype="bar")
fair.save_fig(dirs, "spectratype_cell_type")

# %%
# Step: clones shared between samples (table + heatmap)
# Input: mdata
# Output: tables/18_repertoire_overlap.csv; figures/18_repertoire_overlap_sample.pdf/.png
overlap, _, _ = ir.tl.repertoire_overlap(mdata, "gex:sample", inplace=False)
overlap.to_csv(fair.out(dirs, "tables", "repertoire_overlap.csv"))
grid = ir.pl.repertoire_overlap(mdata, "gex:sample", yticklabels=True, xticklabels=True)
fair.save_fig(dirs, "repertoire_overlap_sample", fig=grid)

# %%
# Step: clonotype clusters by amino acid similarity (TCRdist) and network plot
# Input: mdata
# Output: obs['airr:cc_aa_tcrdist']; figures/18_clonotype_network_tcrdist.pdf/.png
ir.pp.ir_dist(mdata, metric="tcrdist", sequence="aa", cutoff=TCRDIST_CUTOFF)
ir.tl.define_clonotype_clusters(mdata, sequence="aa", metric="tcrdist", receptor_arms="all", dual_ir="any")
ir.tl.clonotype_network(mdata, min_cells=3, sequence="aa", metric="tcrdist")
ir.pl.clonotype_network(mdata, color="gex:cell_type", panel_size=(7, 7))
fair.save_fig(dirs, "clonotype_network_tcrdist")

# %%
# Step: save
# Input: mdata
# Output: objects/18_scirpy.h5mu
mdata.write(fair.out(dirs, "objects", "scirpy.h5mu"))

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
