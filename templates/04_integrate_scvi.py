# 04 scVI / scANVI 整合 + scib-metrics 评估
#
# 【这份代码做什么】
#   用深度生成模型 scVI 从原始 counts 学习去掉批次效应的低维表示（X_scVI），再聚类和画 UMAP。
#   已有部分细胞标签时可以用 scANVI 进一步微调；有标签时还可以用 scib-metrics 定量比较几种整合结果。
#   流程：取高变基因 → 训练 scVI → 看训练曲线 → 潜在空间聚类/UMAP →（可选）scANVI →（可选）scib 评估 → 保存
# 【在流程中的位置】上游 02 → 本步 → 下游 05_annotate（与 03 二选一，05 的 IN_FILE 改成本步输出）
# 【运行环境】sc-scvi（PyTorch；课程 ps10）
# 【输入数据】
#   - results/02_qc_cluster/objects/02_processed.h5ad：layers['counts']（必须是原始 counts）、obs['sample']；scANVI / scib 还需要 obs['cell_type']
# 【输出结果】都在 results/04_scvi/
#   - objects/04_scvi.h5ad：obsm['X_scVI']（可选 X_scANVI）、leiden、X_umap
#   - objects/04_scvi_model/：训练好的模型，可重新加载，不用再训练
#   - figures/：04_training_elbo、04_umap_scvi_sample_leiden、（可选）04_scib_results
#   - tables/：04_cells_per_cluster_sample.csv、（可选）04_scib_metrics.csv
# 【怎么解读】
#   - 训练曲线：ELBO 应平稳下降并收敛；验证集曲线明显回升 = 过拟合，可减少 train() 的 max_epochs
#   - UMAP：按 sample 混合、按细胞类型分开 = 整合良好
#   - scib 表：Batch correction 分数高 = 批次去得干净；Bio conservation 分数高 = 生物差异保留得好；两者要一起看，X_pca 是未整合的基线
# 【注意】
#   - 没有 GPU 也能跑，但细胞多时较慢

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/04_scvi/）
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import scvi
import fair_utils as fair

IN_FILE = "results/02_qc_cluster/objects/02_processed.h5ad"
BATCH_KEY = "sample"
LABEL_KEY = "cell_type"      # 只有跑 scANVI / scib 时需要
N_LATENT = 30
RESOLUTION = 0.8
SEED = 0
scvi.settings.seed = SEED
dirs = fair.setup_step("04_scvi")

# %%
# 步骤：读取数据，只保留高变基因用于建模
# 输入：IN_FILE
# 输出：adata（全部基因）、adata_hvg（高变基因子集）
adata = sc.read_h5ad(IN_FILE)
adata_hvg = adata[:, adata.var.highly_variable].copy()
adata_hvg

# %%
# 步骤：注册数据并训练 scVI（输入必须是原始 counts）
# 输入：adata_hvg.layers['counts']
# 输出：model（训练好的 SCVI）
scvi.model.SCVI.setup_anndata(adata_hvg, layer="counts", batch_key=BATCH_KEY)
model = scvi.model.SCVI(adata_hvg, n_layers=2, n_latent=N_LATENT, gene_likelihood="nb")
model.train()

# %%
# 步骤：训练曲线
# 输入：model.history
# 输出：figures/04_training_elbo.pdf/.png
model.history["elbo_train"].plot()
fair.save_fig(dirs, "training_elbo")

# %%
# 步骤：取潜在空间，聚类和 UMAP
# 输入：model
# 输出：adata.obsm['X_scVI']、obs['leiden']、obsm['X_umap']；figures/04_umap_scvi_sample_leiden.pdf/.png
adata.obsm["X_scVI"] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI", random_state=SEED)
sc.tl.leiden(adata, resolution=RESOLUTION, flavor="igraph", n_iterations=2, directed=False, random_state=SEED)
sc.tl.umap(adata, min_dist=0.3, random_state=SEED)
sc.pl.umap(adata, color=[BATCH_KEY, "leiden"], show=False)
fair.save_fig(dirs, "umap_scvi_sample_leiden")

# %%
# 步骤：每个群里各样本的细胞数
# 输入：adata.obs
# 输出：tables/04_cells_per_cluster_sample.csv
pd.crosstab(adata.obs["leiden"], adata.obs[BATCH_KEY]).to_csv(fair.out(dirs, "tables", "cells_per_cluster_sample.csv"))

# %%
# 可选 scANVI：已有部分或全部细胞类型标签时，在 scVI 基础上微调
# 步骤：从 scVI 初始化 scANVI 并训练
# 输入：model、adata_hvg.obs[LABEL_KEY]（未注释的细胞填 "Unknown"）
# 输出：adata.obsm['X_scANVI']
# adata_hvg.obs[LABEL_KEY] = adata.obs[LABEL_KEY]
# scanvi_model = scvi.model.SCANVI.from_scvi_model(model, adata=adata_hvg, labels_key=LABEL_KEY, unlabeled_category="Unknown")
# scanvi_model.train(max_epochs=20, n_samples_per_label=100)
# adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata_hvg)

# %%
# 可选 scib 评估：需要细胞类型标签
# 步骤：比较 X_pca、X_scVI、X_scANVI 的整合质量
# 输入：adata（obs[BATCH_KEY]、obs[LABEL_KEY]、各 obsm）
# 输出：tables/04_scib_metrics.csv；figures/04_scib_results.pdf/.png
# from scib_metrics.benchmark import Benchmarker
# bm = Benchmarker(adata, batch_key=BATCH_KEY, label_key=LABEL_KEY,
#                  embedding_obsm_keys=["X_pca", "X_scVI", "X_scANVI"], n_jobs=-1)
# bm.benchmark()
# bm.get_results(min_max_scale=False).to_csv(fair.out(dirs, "tables", "scib_metrics.csv"))
# bm.plot_results_table(min_max_scale=False)
# fair.save_fig(dirs, "scib_results")

# %%
# 步骤：保存模型和数据
# 输入：model、adata
# 输出：objects/04_scvi_model/、objects/04_scvi.h5ad
model.save(fair.out(dirs, "objects", "scvi_model"), overwrite=True)
adata.write(fair.out(dirs, "objects", "scvi.h5ad"))

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
