# 01 Read single-cell data and merge samples
#
# [Purpose]
#   Reads the Cell Ranger output of each sample (the three 10x files), merges them into one AnnData object,
#   and adds the sample metadata (condition etc.) to every cell. This is the start of the pipeline; nothing is filtered or normalised here.
#   Workflow: read each sample -> merge (sample suffix added to cell names to avoid duplicates) -> add sample metadata -> save
# [Pipeline position] Start -> downstream 02_qc_cluster
# [Environment] sc-core (course ps6 reading formats, ps7 merging samples)
# [Input]
#   - One 10x folder per sample: matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz (paths in SAMPLES)
#   - data/sample_info.csv: one row per sample; first column sample (matching the names in SAMPLES), other columns are sample metadata; condition is required, batch, sex etc. are optional
# [Output] All in results/01_read_data/
#   - objects/01_merged_raw.h5ad: X = raw counts; obs holds sample and every column of sample_info.csv
#   - tables/01_cells_per_sample.csv: number of cells read per sample
#   - logs/, README.md: run parameters and package versions
# [How to read the results]
#   - cells_per_sample: cell numbers should match the Cell Ranger report (Estimated Number of Cells); a sample with clearly fewer cells needs its sequencing quality checked
#   - The obs summary printed during the run: check each sample's condition is correct; group comparisons in 06 and 07 rely on it
# [Notes]
#   - Files downloaded from GEO often carry a prefix (e.g. GSM123_matrix.mtx.gz); use the prefix= argument or rename them to the standard names
#   - Readers for other formats (h5, h5ad, csv, loom) are in the "Other formats" block; swap it in for the reading block

# %%
# Step: import packages, set the sample table and output folder
# Input: none
# Output: parameter variables; dirs (figures, tables, objects, logs under results/01_read_data/)
from pathlib import Path
import pandas as pd
import scanpy as sc
import anndata as ad
import fair_utils as fair

RAW_DIR = Path("data/raw")
SAMPLE_INFO = "data/sample_info.csv"
SAMPLES = {                      # sample name: folder
    "S1": RAW_DIR / "S1",
    "S2": RAW_DIR / "S2",
}
dirs = fair.setup_step("01_read_data")

# %%
# Step: read each 10x folder
# Input: each folder in SAMPLES
# Output: adatas (dict, sample name -> AnnData)
adatas = {}
for name, folder in SAMPLES.items():
    adatas[name] = sc.read_10x_mtx(folder, var_names="gene_symbols", cache=True)
    adatas[name].var_names_make_unique()
adatas

# %%
# Other formats (replace the block above as needed):
# sc.read_10x_h5("sample_filtered_feature_bc_matrix.h5")
# sc.read_h5ad("sample.h5ad")
# sc.read_csv("counts.csv").T      # transpose when the file is genes x cells to get cells x genes
# sc.read_text("counts.txt").T
# sc.read_10x_mtx(folder, prefix="GSM123_")

# %%
# Step: merge all samples, adding the sample name to cell names to avoid duplicates
# Input: adatas
# Output: adata (cells x genes, obs['sample'])
adata = ad.concat(adatas, label="sample", index_unique="-")
adata

# %%
# Step: add the sample metadata (condition etc.) to every cell
# Input: adata.obs['sample'], SAMPLE_INFO
# Output: adata.obs gains the columns of sample_info
sample_info = pd.read_csv(SAMPLE_INFO, index_col="sample")
adata.obs = adata.obs.join(sample_info, on="sample")
adata.obs.groupby("sample", observed=True).first()

# %%
# Step: save the number of cells per sample and the merged raw data
# Input: adata
# Output: tables/01_cells_per_sample.csv, objects/01_merged_raw.h5ad
adata.obs["sample"].value_counts().to_csv(fair.out(dirs, "tables", "cells_per_sample.csv"))
adata.write(fair.out(dirs, "objects", "merged_raw.h5ad"))

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
