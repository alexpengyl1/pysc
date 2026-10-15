# 08c 找细胞类型特异的 regulon 并作图
#
# 【这份代码做什么】
#   把 08b 算出的 regulon 活性矩阵读回 AnnData，按细胞类型找出各自特异活跃的转录因子，并画矩阵图和聚类热图。
# 【在流程中的位置】上游 08b → 本步（终点分析）
# 【运行环境】sc-core（课程 ps13；课程只画了 clustermap，这里增加了按细胞类型的汇总）
# 【输入数据】
#   - results/08_scenic/tables/08_auc_matrix.csv（08b 输出）
#   - results/08_scenic/tables/08_cell_info.csv（08a 输出）
# 【输出结果】都在 results/08_scenic/
#   - tables/08_regulon_markers_cell_type.csv：每个细胞类型特异的 regulon（Wilcoxon，一对其余）
#   - figures/：08_regulon_matrixplot_cell_type、08_auc_clustermap
# 【怎么解读】
#   - matrixplot：每行一个细胞类型、每列一个 regulon，颜色 = 平均活性（按列标准化）；某行特别亮的 regulon = 该细胞类型的身份或状态调控因子
#   - regulon_markers 表：scores 高、pvals_adj 小的 regulon 最特异
#   - clustermap：细胞按 regulon 活性聚类，能否按细胞类型分开反映调控程序的差异

# %%
# 步骤：导入包，设置路径和输出文件夹
# 输入：无
# 输出：dirs（results/08_scenic/）
import pandas as pd
import seaborn as sns
import scanpy as sc
import anndata as ad
import fair_utils as fair

AUC_FILE = "results/08_scenic/tables/08_auc_matrix.csv"
CELL_INFO_FILE = "results/08_scenic/tables/08_cell_info.csv"
N_TOP = 3
dirs = fair.setup_step("08_scenic")

# %%
# 步骤：读取 AUC 矩阵，做成以 regulon 为特征的 AnnData
# 输入：AUC_FILE、CELL_INFO_FILE
# 输出：auc_adata（细胞 x regulon，obs['cell_type']）
auc = pd.read_csv(AUC_FILE, index_col=0)
cell_info = pd.read_csv(CELL_INFO_FILE, index_col=0)
auc_adata = ad.AnnData(auc)
auc_adata.obs["cell_type"] = cell_info.loc[auc.index, "cell_type"].astype("category")
auc_adata

# %%
# 步骤：每个细胞类型特异的 regulon（Wilcoxon，一对其余）
# 输入：auc_adata
# 输出：tables/08_regulon_markers_cell_type.csv
sc.tl.rank_genes_groups(auc_adata, "cell_type", method="wilcoxon")
sc.get.rank_genes_groups_df(auc_adata, group=None).to_csv(
    fair.out(dirs, "tables", "regulon_markers_cell_type.csv"), index=False)

# %%
# 步骤：各细胞类型前 N_TOP 个特异 regulon 的活性矩阵图
# 输入：auc_adata
# 输出：figures/08_regulon_matrixplot_cell_type.pdf/.png
sc.pl.rank_genes_groups_matrixplot(auc_adata, n_genes=N_TOP, standard_scale="var", dendrogram=False, show=False)
fair.save_fig(dirs, "regulon_matrixplot_cell_type")

# %%
# 步骤：全部细胞 x regulon 的聚类热图（课程 ps13 的图）
# 输入：auc
# 输出：figures/08_auc_clustermap.pdf/.png
grid = sns.clustermap(auc, figsize=(12, 12), yticklabels=False)
fair.save_fig(dirs, "auc_clustermap", fig=grid)

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info_08c.json、logs/pip_freeze_08c.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()}, script="08c")
