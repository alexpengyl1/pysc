# 生成模板测试用的小数据集（由 run_all_tests.py 在 sc-core 环境里自动调用，也可单独运行）
# 环境：sc-core（需要联网下载 PBMC 3k）
# 输入：scanpy 自带的公开数据 pbmc3k（原始 counts）和 pbmc3k_processed（官方细胞类型注释）
# 输出（tests/test_project/data/）：
#   raw/S1–S4/           4 个样本的 10x 三件套（S1、S2 = ctrl，S3、S4 = stim）
#   sample_info.csv      样本 -> condition
#   test_cell_types.csv  barcode -> cell_type（05 注释步骤用它代替人工对照表）
#   geneset.gs           scDRS 测试基因集（两个假“性状”：单核细胞基因、NK 基因）
# 人为设计的信号（只为让下游步骤有结果，不代表真实生物学）：
#   stim 样本中干扰素基因 x4（上调）、HLA/S100 基因保留 1/4（下调）

# %%
# 步骤：导入包，设置参数
# 输入：无
# 输出：参数变量
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
CELL_TYPE_NAMES = {                       # 官方注释 -> 文件名友好的名字
    "CD4 T cells": "CD4_T", "CD14+ Monocytes": "CD14_Mono", "B cells": "B", "CD8 T cells": "CD8_T",
    "NK cells": "NK", "FCGR3A+ Monocytes": "FCGR3A_Mono", "Dendritic cells": "DC", "Megakaryocytes": "Mk",
}
SEED = 0
rng = np.random.default_rng(SEED)

# %%
# 步骤：下载原始 counts 和官方注释，只保留有注释的细胞
# 输入：网络（10x Genomics、cellxgene 示例数据）
# 输出：adata（细胞 x 基因，原始 counts，obs['cell_type']）
adata = sc.datasets.pbmc3k()
annotated = sc.datasets.pbmc3k_processed()
adata = adata[annotated.obs_names].copy()
adata.obs["cell_type"] = annotated.obs["louvain"].map(CELL_TYPE_NAMES).values
adata.var_names_make_unique()
adata

# %%
# 步骤：随机把细胞分到 4 个样本
# 输入：adata
# 输出：adata.obs['sample']、adata.obs['condition']
adata.obs["sample"] = rng.choice(list(SAMPLES), size=adata.n_obs)
adata.obs["condition"] = adata.obs["sample"].map(SAMPLES)
adata.obs["sample"].value_counts()

# %%
# 步骤：在 stim 样本里上调干扰素基因、下调 HLA/S100 基因
# 输入：adata
# 输出：adata.X（修改后的 counts）
counts = adata.X.toarray()
stim = (adata.obs["condition"] == "stim").values
up = adata.var_names.isin(UP_GENES)
down = adata.var_names.isin(DOWN_GENES)
counts[np.ix_(stim, up)] = rng.poisson(counts[np.ix_(stim, up)] * 4 + 0.5)
counts[np.ix_(stim, down)] = rng.binomial(counts[np.ix_(stim, down)].astype(int), 0.25)
adata.X = sp.csr_matrix(counts)
print("上调基因：", list(adata.var_names[up]), "下调基因：", list(adata.var_names[down]))

# %%
# 步骤：每个样本写成 10x 三件套（matrix.mtx.gz / features.tsv.gz / barcodes.tsv.gz）
# 输入：adata
# 输出：data/raw/S1–S4/
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
# 步骤：样本分组表和细胞类型表
# 输入：SAMPLES、adata.obs
# 输出：data/sample_info.csv、data/test_cell_types.csv
pd.DataFrame({"sample": list(SAMPLES), "condition": list(SAMPLES.values())}).to_csv(DATA_DIR / "sample_info.csv", index=False)
adata.obs[["cell_type"]].rename_axis("barcode").to_csv(DATA_DIR / "test_cell_types.csv")

# %%
# 步骤：scDRS 测试基因集（格式：TRAIT<TAB>GENE1:权重,GENE2:权重）
# 输入：无
# 输出：data/geneset.gs
gene_sets = {
    "MONO_TEST": ["CD14", "LYZ", "S100A12", "FCN1", "VCAN", "CSF1R", "CD68", "FCGR3A", "TYROBP", "CST3"],
    "NK_TEST": ["NKG7", "GNLY", "GZMB", "PRF1", "KLRD1", "KLRF1", "FCGR3A", "CTSW", "GZMA", "CST7"],
}
gs = pd.DataFrame({"TRAIT": list(gene_sets), "GENESET": [",".join(f"{g}:1" for g in genes) for genes in gene_sets.values()]})
gs.to_csv(DATA_DIR / "geneset.gs", sep="\t", index=False)
print("测试数据已写入", DATA_DIR.resolve())
