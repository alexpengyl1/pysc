# 03 多样本整合：Harmony（首选）/ BBKNN / ingest
#
# 【这份代码做什么】
#   去除样本之间的批次效应，让同一种细胞不管来自哪个样本都聚在一起，然后重新聚类和画 UMAP。
#   提供三种方法，只选一种运行：
#     A Harmony：在 PCA 空间里校正批次（最常用）
#     B BBKNN：构建批次平衡的邻居图，不改 PCA
#     C ingest：有已经注释好的参考数据时，把参考的标签和 UMAP 投到新数据上
# 【在流程中的位置】上游 02 → 本步 → 下游 05_annotate（与 04 scVI 二选一）
# 【运行环境】sc-core（课程 ps7）
# 【输入数据】
#   - results/02_qc_cluster/objects/02_processed.h5ad（需要 obsm['X_pca']、obs['sample']）
#   - 方法 C 还需要一份已注释的参考 h5ad
# 【输出结果】都在 results/03_integrate/
#   - objects/03_integrated.h5ad：obsm['X_pca_harmony']（方法 A），重算的 leiden 和 X_umap
#   - figures/：03_umap_before_integration、03_umap_after_harmony（或 bbknn / ingest）、03_markers_dotplot_leiden
#   - tables/：03_markers_leiden.csv、03_cells_per_cluster_sample.csv（每个群里各样本的细胞数）
# 【怎么解读】
#   - 整合前后 UMAP 对比：整合后按 sample 上色应混在一起，按 marker 上色仍然分开
#   - 同一细胞类型仍按样本拆开 = 校正不足；不同细胞类型被糊在一起 = 校正过度
#   - cells_per_cluster_sample 表：某个群只来自一个样本，要判断是批次残留还是真实的样本特异细胞群（例如只在处理组出现的状态）

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/03_integrate/）
import pandas as pd
import scanpy as sc
import fair_utils as fair

IN_FILE = "results/02_qc_cluster/objects/02_processed.h5ad"
METHOD = "harmony"       # "harmony" / "bbknn" / "ingest"，只是记录用，实际运行哪段由你选择
BATCH_KEY = "sample"
N_PCS = 30
RESOLUTION = 0.8
dirs = fair.setup_step("03_integrate")

# %%
# 步骤：读取上一步结果，画整合前的 UMAP
# 输入：IN_FILE
# 输出：adata；figures/03_umap_before_integration.pdf/.png
adata = sc.read_h5ad(IN_FILE)
sc.pl.umap(adata, color=[BATCH_KEY, "leiden"], show=False)
fair.save_fig(dirs, "umap_before_integration")

# %%
# 方法 A — Harmony
# 步骤：在 PCA 空间校正批次
# 输入：adata.obsm['X_pca']
# 输出：adata.obsm['X_pca_harmony']
sc.external.pp.harmony_integrate(adata, key=BATCH_KEY, basis="X_pca", adjusted_basis="X_pca_harmony")

# %%
# 步骤：用校正后的空间重算邻居图、聚类、UMAP
# 输入：adata.obsm['X_pca_harmony']
# 输出：adata.obs['leiden']、adata.obsm['X_umap']；figures/03_umap_after_harmony.pdf/.png
sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_pcs=N_PCS, random_state=0)
sc.tl.leiden(adata, resolution=RESOLUTION, flavor="igraph", n_iterations=2, directed=False, random_state=0)
sc.tl.umap(adata, random_state=0)
sc.pl.umap(adata, color=[BATCH_KEY, "leiden"], show=False)
fair.save_fig(dirs, "umap_after_harmony")

# %%
# 方法 B — BBKNN（替代 A：直接构建批次平衡的邻居图，不改表达矩阵和 PCA）
# 步骤：BBKNN 邻居图 + 聚类 + UMAP
# 输入：adata.obsm['X_pca']
# 输出：adata.obs['leiden']、adata.obsm['X_umap']；figures/03_umap_after_bbknn.pdf/.png
# sc.external.pp.bbknn(adata, batch_key=BATCH_KEY, n_pcs=N_PCS)
# sc.tl.leiden(adata, resolution=RESOLUTION, flavor="igraph", n_iterations=2, directed=False, random_state=0)
# sc.tl.umap(adata, random_state=0)
# sc.pl.umap(adata, color=[BATCH_KEY, "leiden"], show=False)
# fair.save_fig(dirs, "umap_after_bbknn")

# %%
# 方法 C — ingest（有注释好的参考数据时，把参考的标签和 UMAP 投到新数据上）
# 步骤：两份数据取共同基因，参考数据算 PCA/邻居/UMAP，再映射
# 输入：adata_ref（已注释，obs['cell_type']），adata
# 输出：adata.obs['cell_type']（来自参考）、adata.obsm['X_umap']；figures/03_umap_ingest_cell_type.pdf/.png
# adata_ref = sc.read_h5ad("data/reference.h5ad")
# genes = adata_ref.var_names.intersection(adata.var_names)
# adata_ref = adata_ref[:, genes].copy()
# adata = adata[:, genes].copy()
# sc.pp.pca(adata_ref)
# sc.pp.neighbors(adata_ref)
# sc.tl.umap(adata_ref)
# sc.tl.ingest(adata, adata_ref, obs="cell_type")
# sc.pl.umap(adata, color="cell_type", show=False)
# fair.save_fig(dirs, "umap_ingest_cell_type")

# %%
# 步骤：每个群里各样本的细胞数（检查批次残留）
# 输入：adata.obs
# 输出：tables/03_cells_per_cluster_sample.csv
cluster_sample = pd.crosstab(adata.obs["leiden"], adata.obs[BATCH_KEY])
cluster_sample.to_csv(fair.out(dirs, "tables", "cells_per_cluster_sample.csv"))
cluster_sample

# %%
# 步骤：整合后的 marker 基因
# 输入：adata
# 输出：tables/03_markers_leiden.csv；figures/03_markers_dotplot_leiden.pdf/.png
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=False)
sc.get.rank_genes_groups_df(adata, group=None).to_csv(fair.out(dirs, "tables", "markers_leiden.csv"), index=False)
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, dendrogram=False, show=False)
fair.save_fig(dirs, "markers_dotplot_leiden")

# %%
# 步骤：保存
# 输入：adata
# 输出：objects/03_integrated.h5ad
adata.write(fair.out(dirs, "objects", "integrated.h5ad"))

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
