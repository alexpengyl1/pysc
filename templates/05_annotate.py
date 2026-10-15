# 05 细胞类型注释：marker 手动注释 + CellTypist 自动注释交叉核对
#
# 【这份代码做什么】
#   给每个聚类群命名细胞类型。先看已知 marker 基因在各群的表达，再用 CellTypist 预训练模型自动注释，
#   两者对照后手动写出 群 → 细胞类型 的对照表，生成后续所有分析都要用的 obs['cell_type']。
#   流程：marker dotplot → CellTypist 自动注释 → 与聚类对照 → 手动填写对照表 → 注释结果作图 → 保存
# 【在流程中的位置】上游 03 或 04 → 本步 → 下游 06–24 全部分析（都读 05_annotated.h5ad）
# 【运行环境】sc-core（课程 ps2/ps6 手动注释、ps9 CellTypist）
# 【输入数据】
#   - results/03_integrate/objects/03_integrated.h5ad（或 02 / 04 的结果）；X = log1p（每细胞总量 1e4，全部基因）
#   - 顶部 MARKERS：你的组织里各细胞类型的已知 marker（数据里没有的基因会自动去掉并打印出来）
#   - CellTypist 模型：第一次运行会联网下载
# 【输出结果】都在 results/05_annotate/
#   - objects/05_annotated.h5ad：obs['cell_type']（最终注释）、obs['celltypist_*']（自动注释和置信度）
#   - figures/：05_markers_dotplot_leiden、05_umap_celltypist、05_celltypist_vs_leiden_dotplot、05_umap_cell_type、05_markers_dotplot_cell_type
#   - tables/：05_cluster_annotation.csv（群 → 细胞类型对照表）、05_celltypist_vs_leiden.csv、05_cells_per_type_sample.csv
# 【怎么解读】
#   - marker dotplot：某群里一组 marker 同时高表达（点大且深）= 该细胞类型
#   - CellTypist：majority_voting（按群多数投票）比逐细胞 predicted_labels 更稳；conf_score 低的细胞要人工复查
#   - celltypist_vs_leiden：一个群对应多个 CellTypist 标签，可能是群分得太粗；两种结果矛盾的群，回到 marker 判断
#   - cells_per_type_sample：各样本的细胞组成，是 06 组成分析的原始数据
# 【注意】
#   - 对照表要看完前面的图后手动填写，不能直接用默认值
#   - CellTypist 模型大多是人类的；小鼠肠道可用 Adult_Mouse_Gut，其他小鼠数据要谨慎

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/05_annotate/）
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models
import fair_utils as fair

IN_FILE = "results/03_integrate/objects/03_integrated.h5ad"
CLUSTER_KEY = "leiden"
CELLTYPIST_MODEL = "Immune_All_Low.pkl"    # 其他：Immune_All_High、Cells_Intestinal_Tract、Adult_Mouse_Gut …
MARKERS = {                                # 按自己的组织改
    "T": ["CD3D", "CD3E"],
    "B": ["MS4A1", "CD79A"],
    "Myeloid": ["CD14", "LYZ"],
    "Epithelial": ["EPCAM", "KRT8"],
    "Fibroblast": ["COL1A1", "DCN"],
    "Endothelial": ["PECAM1", "VWF"],
}
dirs = fair.setup_step("05_annotate")

# %%
# 步骤：读取数据，只保留数据里存在的 marker 基因
# 输入：IN_FILE、MARKERS
# 输出：adata、markers（去掉不存在基因后的字典）；打印被去掉的基因
adata = sc.read_h5ad(IN_FILE)
markers = {cell: [g for g in genes if g in adata.var_names] for cell, genes in MARKERS.items()}
markers = {cell: genes for cell, genes in markers.items() if genes}
print("数据中没有的 marker：", [g for genes in MARKERS.values() for g in genes if g not in adata.var_names])

# %%
# 步骤：各 leiden 群的 marker dotplot
# 输入：adata、markers
# 输出：figures/05_markers_dotplot_leiden.pdf/.png
sc.pl.dotplot(adata, markers, groupby=CLUSTER_KEY, show=False)
fair.save_fig(dirs, "markers_dotplot_leiden")

# %%
# 步骤：下载 CellTypist 模型（第一次运行需要），查看可用模型
# 输入：无
# 输出：本地模型缓存；打印模型说明
models.download_models(force_update=False)
models.models_description()

# %%
# 步骤：CellTypist 注释（多数投票）
# 输入：adata（X = log1p）
# 输出：predictions；adata.obs['celltypist_label']、['celltypist_majority']、['celltypist_conf']
predictions = celltypist.annotate(adata, model=CELLTYPIST_MODEL, majority_voting=True)
adata.obs["celltypist_label"] = predictions.predicted_labels["predicted_labels"]
adata.obs["celltypist_majority"] = predictions.predicted_labels["majority_voting"]
adata.obs["celltypist_conf"] = predictions.probability_matrix.max(axis=1)

# %%
# 步骤：CellTypist 结果画在 UMAP 上
# 输入：adata.obs
# 输出：figures/05_umap_celltypist.pdf/.png
sc.pl.umap(adata, color=[CLUSTER_KEY, "celltypist_majority", "celltypist_conf"], legend_loc="on data", show=False)
fair.save_fig(dirs, "umap_celltypist")

# %%
# 步骤：CellTypist 结果和聚类的对照（图 + 交叉表）
# 输入：predictions、adata.obs
# 输出：figures/05_celltypist_vs_leiden_dotplot.pdf/.png；tables/05_celltypist_vs_leiden.csv
celltypist.dotplot(predictions, use_as_reference=CLUSTER_KEY, use_as_prediction="majority_voting", show=False)
fair.save_fig(dirs, "celltypist_vs_leiden_dotplot")
pd.crosstab(adata.obs[CLUSTER_KEY], adata.obs["celltypist_majority"]).to_csv(fair.out(dirs, "tables", "celltypist_vs_leiden.csv"))

# %%
# 步骤：手动写 cluster -> 细胞类型 对照表（看完上面的图后填写），并保存对照表
# 输入：adata.obs[CLUSTER_KEY]
# 输出：adata.obs['cell_type']；tables/05_cluster_annotation.csv
annotation = {
    "0": "T",
    "1": "B",
    # ...
}
adata.obs["cell_type"] = adata.obs[CLUSTER_KEY].map(annotation).astype("category")
pd.Series(annotation, name="cell_type").rename_axis(CLUSTER_KEY).to_csv(fair.out(dirs, "tables", "cluster_annotation.csv"))

# %%
# 步骤：注释结果 UMAP 和 marker dotplot
# 输入：adata
# 输出：figures/05_umap_cell_type.pdf/.png、05_markers_dotplot_cell_type.pdf/.png
sc.pl.umap(adata, color="cell_type", legend_loc="on data", show=False)
fair.save_fig(dirs, "umap_cell_type")
sc.pl.dotplot(adata, markers, groupby="cell_type", show=False)
fair.save_fig(dirs, "markers_dotplot_cell_type")

# %%
# 步骤：每个样本各细胞类型的细胞数
# 输入：adata.obs
# 输出：tables/05_cells_per_type_sample.csv
pd.crosstab(adata.obs["cell_type"], adata.obs["sample"]).to_csv(fair.out(dirs, "tables", "cells_per_type_sample.csv"))

# %%
# 步骤：保存
# 输入：adata
# 输出：objects/05_annotated.h5ad
adata.write(fair.out(dirs, "objects", "annotated.h5ad"))

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
