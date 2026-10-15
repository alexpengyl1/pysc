# 01 读取单细胞数据并合并样本
#
# 【这份代码做什么】
#   把每个样本的 Cell Ranger 输出（10x 三件套）读进来，合并成一个 AnnData 对象，
#   并把样本分组信息（condition 等）写到每个细胞上。这是整个流程的起点，不做任何过滤和标准化。
#   流程：逐个读取样本 → 合并（细胞名加样本后缀防重复）→ 合并样本元数据 → 保存
# 【在流程中的位置】起点 → 下游 02_qc_cluster
# 【运行环境】sc-core（课程 ps6 多种格式读取、ps7 多样本合并）
# 【输入数据】
#   - 每个样本一个 10x 文件夹：matrix.mtx.gz、features.tsv.gz、barcodes.tsv.gz（路径写在 SAMPLES）
#   - data/sample_info.csv：每个样本一行，第一列 sample（与 SAMPLES 的名字一致），其余列是样本信息，至少要有 condition，可加 batch、sex 等
# 【输出结果】都在 results/01_read_data/
#   - objects/01_merged_raw.h5ad：X = 原始 counts；obs 里有 sample 和 sample_info.csv 的全部列
#   - tables/01_cells_per_sample.csv：每个样本读入的细胞数
#   - logs/、README.md：运行参数和包版本
# 【怎么解读】
#   - 看 cells_per_sample：每个样本的细胞数应与 Cell Ranger 报告（Estimated Number of Cells）一致；某个样本明显偏少要回查测序质量
#   - 运行时打印的 obs 汇总表：确认每个样本的 condition 写对了，后面 06、07 的组间比较都靠它
# 【注意】
#   - GEO 下载的文件常带前缀（如 GSM123_matrix.mtx.gz），用 prefix= 参数或改名为标准三件套
#   - 其他格式（h5、h5ad、csv、loom）的读法在「其他格式」块里，替换读取块即可

# %%
# 步骤：导入包，设置样本表和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/01_read_data/ 下的 figures、tables、objects、logs）
from pathlib import Path
import pandas as pd
import scanpy as sc
import anndata as ad
import fair_utils as fair

RAW_DIR = Path("data/raw")
SAMPLE_INFO = "data/sample_info.csv"
SAMPLES = {                      # 样本名: 文件夹
    "S1": RAW_DIR / "S1",
    "S2": RAW_DIR / "S2",
}
dirs = fair.setup_step("01_read_data")

# %%
# 步骤：逐个读取 10x 文件夹
# 输入：SAMPLES 里的每个文件夹
# 输出：adatas（字典，样本名 -> AnnData）
adatas = {}
for name, folder in SAMPLES.items():
    adatas[name] = sc.read_10x_mtx(folder, var_names="gene_symbols", cache=True)
    adatas[name].var_names_make_unique()
adatas

# %%
# 其他格式（按需替换上一块）：
# sc.read_10x_h5("sample_filtered_feature_bc_matrix.h5")
# sc.read_h5ad("sample.h5ad")
# sc.read_csv("counts.csv").T      # 文件是 基因 x 细胞 时要转置成 细胞 x 基因
# sc.read_text("counts.txt").T
# sc.read_10x_mtx(folder, prefix="GSM123_")

# %%
# 步骤：合并所有样本，细胞名加上样本后缀避免重复
# 输入：adatas
# 输出：adata（细胞 x 基因，obs['sample']）
adata = ad.concat(adatas, label="sample", index_unique="-")
adata

# %%
# 步骤：把样本元数据（condition 等）合并到每个细胞
# 输入：adata.obs['sample']、SAMPLE_INFO
# 输出：adata.obs 新增 sample_info 里的列
sample_info = pd.read_csv(SAMPLE_INFO, index_col="sample")
adata.obs = adata.obs.join(sample_info, on="sample")
adata.obs.groupby("sample", observed=True).first()

# %%
# 步骤：保存每个样本的细胞数和合并后的原始数据
# 输入：adata
# 输出：tables/01_cells_per_sample.csv、objects/01_merged_raw.h5ad
adata.obs["sample"].value_counts().to_csv(fair.out(dirs, "tables", "cells_per_sample.csv"))
adata.write(fair.out(dirs, "objects", "merged_raw.h5ad"))

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
