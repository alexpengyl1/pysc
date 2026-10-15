# 12 PHATE 降维可视化
#
# 【这份代码做什么】
#   用 PHATE 把细胞画到二维平面上。和 UMAP 相比，PHATE 更能保留连续的分化或时间序列结构（分支、过渡态）。
#   只做可视化，不产生拟时序或统计结果。
# 【在流程中的位置】上游 05 → 本步（终点分析）
# 【运行环境】sc-core（课程 ps17）
# 【输入数据】
#   - results/05_annotate/objects/05_annotated.h5ad（X = log1p）
# 【输出结果】都在 results/12_phate/
#   - tables/12_phate_coordinates.csv：每个细胞的 PHATE1 / PHATE2 坐标
#   - figures/12_phate_cell_type_sample
# 【怎么解读】
#   - 按细胞类型 / 时间点 / 样本上色，看是否有连续的轨迹和分支；分支端点通常是终末细胞类型
#   - 参数：KNN 小 → 局部结构更细；大数据集（> 10 万细胞）可把 KNN 调到 30–40；T 控制平滑程度，越大越平滑
# 【注意】
#   - 图上的距离和形状只用于直观展示，结论要配合拟时序（11）或 RNA velocity（13）

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/12_phate/）
import pandas as pd
import scanpy as sc
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
KNN = 5
T = "auto"
COLOR_BY = ["cell_type", "sample"]
SEED = 0
dirs = fair.setup_step("12_phate")

# %%
# 步骤：读取数据
# 输入：IN_FILE
# 输出：adata
adata = sc.read_h5ad(IN_FILE)

# %%
# 步骤：计算 PHATE 嵌入并保存坐标
# 输入：adata.X
# 输出：adata.obsm['X_phate']；tables/12_phate_coordinates.csv
sc.external.tl.phate(adata, k=KNN, t=T, n_components=2, n_jobs=-2, random_state=SEED)
pd.DataFrame(adata.obsm["X_phate"], index=adata.obs_names, columns=["PHATE1", "PHATE2"]).to_csv(
    fair.out(dirs, "tables", "phate_coordinates.csv"))

# %%
# 步骤：画 PHATE 图
# 输入：adata.obsm['X_phate']
# 输出：figures/12_phate_cell_type_sample.pdf/.png
sc.external.pl.phate(adata, color=COLOR_BY, show=False)
fair.save_fig(dirs, "phate_cell_type_sample")

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
