# Build a small dataset for testing the templates (called automatically by run_all_tests.py in sc-core, or run on its own)
# Environment: sc-core (needs internet to download PBMC 3k)
# Input: scanpy's public datasets pbmc3k (raw counts) and pbmc3k_processed (official cell type annotation)
# Output (tests/test_project/data/):
#   raw/S1-S4/           10x files for 4 samples (S1, S2 = ctrl; S3, S4 = stim)
#   sample_info.csv      sample -> condition
#   test_cell_types.csv  barcode -> cell_type (used by the 05 annotation step instead of a hand-written mapping)
#   geneset.gs           scDRS test gene sets (two fake "traits": monocyte genes, NK genes)
# Artificial signal (only so downstream steps have something to find, not real biology):
#   in stim samples interferon genes x4 (up), HLA/S100 genes kept at 1/4 (down)

# %%
# Step: import packages, set parameters
# Input: none
# Output: parameter variables
import gzip
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import scanpy as sc

DATA_DIR = Path("data")
SAMPLES = {"S1": "ctrl", "S2": "ctrl", "S3": "stim", "S4": "stim"}
UP_GENES = ["ISG15", "IFI6", "IFIT3", "IFITM3", "MX1", "OAS1", "ISG20", "IFI44L"]
DOWN_GENES = ["HLA-DRA", "HLA-DRB1", "CD74", "S100A8", "S100A9"]
CELL_TYPE_NAMES = {                       # official label -> file-name-friendly name
    "CD4 T cells": "CD4_T", "CD14+ Monocytes": "CD14_Mono", "B cells": "B", "CD8 T cells": "CD8_T",
    "NK cells": "NK", "FCGR3A+ Monocytes": "FCGR3A_Mono", "Dendritic cells": "DC", "Megakaryocytes": "Mk",
}
SEED = 0
rng = np.random.default_rng(SEED)

# %%
# Step: download raw counts and the official annotation, keep annotated cells only
# Input: internet (10x Genomics, cellxgene example data)
# Output: adata (cells x genes, raw counts, obs['cell_type'])
adata = sc.datasets.pbmc3k()
annotated = sc.datasets.pbmc3k_processed()
adata = adata[annotated.obs_names].copy()
adata.obs["cell_type"] = annotated.obs["louvain"].map(CELL_TYPE_NAMES).values
adata.var_names_make_unique()
adata

# %%
# Step: randomly assign cells to 4 samples
# Input: adata
# Output: adata.obs['sample'], adata.obs['condition']
adata.obs["sample"] = rng.choice(list(SAMPLES), size=adata.n_obs)
adata.obs["condition"] = adata.obs["sample"].map(SAMPLES)
adata.obs["sample"].value_counts()

# %%
# Step: up-regulate interferon genes and down-regulate HLA/S100 genes in stim samples
# Input: adata
# Output: adata.X (modified counts)
counts = adata.X.toarray()
stim = (adata.obs["condition"] == "stim").values
up = adata.var_names.isin(UP_GENES)
down = adata.var_names.isin(DOWN_GENES)
counts[np.ix_(stim, up)] = rng.poisson(counts[np.ix_(stim, up)] * 4 + 0.5)
counts[np.ix_(stim, down)] = rng.binomial(counts[np.ix_(stim, down)].astype(int), 0.25)
adata.X = sp.csr_matrix(counts)
print("Up-regulated genes:", list(adata.var_names[up]), "Down-regulated genes:", list(adata.var_names[down]))

# %%
# Step: write each sample as 10x files (matrix.mtx.gz / features.tsv.gz / barcodes.tsv.gz)
# Input: adata
# Output: data/raw/S1-S4/
for sample in SAMPLES:
    sub = adata[adata.obs["sample"] == sample]
    out_dir = DATA_DIR / "raw" / sample
    out_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_dir / "matrix.mtx.gz", "wb") as f:
        scipy.io.mmwrite(f, sp.csr_matrix(sub.X.T).astype(np.int32))
    pd.DataFrame({"id": sub.var["gene_ids"], "name": sub.var_names, "type": "Gene Expression"}).to_csv(
        out_dir / "features.tsv.gz", sep="\t", header=False, index=False)
    pd.Series(sub.obs_names).to_csv(out_dir / "barcodes.tsv.gz", header=False, index=False)

# %%
# Step: sample table and cell type table
# Input: SAMPLES, adata.obs
# Output: data/sample_info.csv, data/test_cell_types.csv
pd.DataFrame({"sample": list(SAMPLES), "condition": list(SAMPLES.values())}).to_csv(DATA_DIR / "sample_info.csv", index=False)
adata.obs[["cell_type"]].rename_axis("barcode").to_csv(DATA_DIR / "test_cell_types.csv")

# %%
# Step: scDRS test gene sets (format: TRAIT<TAB>GENE1:weight,GENE2:weight)
# Input: none
# Output: data/geneset.gs
gene_sets = {
    "MONO_TEST": ["CD14", "LYZ", "S100A12", "FCN1", "VCAN", "CSF1R", "CD68", "FCGR3A", "TYROBP", "CST3"],
    "NK_TEST": ["NKG7", "GNLY", "GZMB", "PRF1", "KLRD1", "KLRF1", "FCGR3A", "CTSW", "GZMA", "CST7"],
}
gs = pd.DataFrame({"TRAIT": list(gene_sets), "GENESET": [",".join(f"{g}:1" for g in genes) for genes in gene_sets.values()]})
gs.to_csv(DATA_DIR / "geneset.gs", sep="\t", index=False)
print("Test data written to", DATA_DIR.resolve())
