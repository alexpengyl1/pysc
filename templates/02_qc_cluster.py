# 02 质控、标准化、高变基因、PCA、聚类、UMAP、marker 基因
#
# 【这份代码做什么】
#   对合并后的原始数据做 scanpy 标准预处理：去掉低质量细胞和双细胞，标准化，降维聚类，
#   并找出每个群的 marker 基因，为后面的整合和注释做准备。
#   流程：QC 指标 → 看图定阈值 → 过滤 → Scrublet 去双细胞 → 标准化 + log1p → 高变基因 → PCA → 邻居图 → Leiden 聚类 → UMAP → marker
# 【在流程中的位置】上游 01 → 本步 → 下游 03（Harmony/BBKNN）或 04（scVI）
# 【运行环境】sc-core（课程 ps2 scanpy 标准流程、ps7）
# 【输入数据】
#   - results/01_read_data/objects/01_merged_raw.h5ad（X = 原始 counts，obs['sample']）
# 【输出结果】都在 results/02_qc_cluster/
#   - objects/02_processed.h5ad：X = log1p 标准化（全部基因，未 scale）；layers['counts'] = 原始 counts；layers['log1p'] = 同 X；
#     obs 有 QC 指标、predicted_doublet、leiden；obsm 有 X_pca、X_umap
#   - figures/：02_qc_violin_by_sample、02_qc_scatter_counts_genes、02_hvg、02_pca_variance_ratio、02_umap_leiden_sample、02_markers_dotplot_leiden（各有 .pdf 和 .png）
#   - tables/：02_qc_cells_kept.csv（每个样本过滤前后细胞数）、02_doublets_per_sample.csv、02_markers_leiden.csv
# 【怎么解读】
#   - QC 小提琴图：n_genes_by_counts 过低 = 空液滴或碎片，过高 = 可能是双细胞；pct_counts_mt 高 = 细胞破损。阈值看完图再改顶部参数
#   - qc_cells_kept 表：某个样本被过滤掉的比例远高于其他样本，说明该样本质量差，需要单独检查
#   - 双细胞比例一般随上样细胞数增加，约 1–8%
#   - PCA 方差图：曲线变平的拐点之后 PC 贡献很小，N_PCS 取拐点附近
#   - UMAP 按 sample 上色：如果样本各自成团，说明有批次效应，下一步要整合
#   - marker dotplot：点大小 = 群内表达该基因的细胞比例，颜色 = 平均表达；用来判断每个 leiden 群是什么细胞
# 【注意】
#   - 阈值（MIN_GENES、MAX_PCT_MT 等）因组织而异，不要照抄默认值；小鼠线粒体前缀用 mt-

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/02_qc_cluster/）
import pandas as pd
import scanpy as sc
import fair_utils as fair

IN_FILE = "results/01_read_data/objects/01_merged_raw.h5ad"
MT_PREFIX = "MT-"        # 人 "MT-"，小鼠 "mt-"
MIN_GENES = 200          # 课程 PBMC：200–2500 基因、线粒体 <5%；肿瘤组织：>500 基因、线粒体 <20%
MAX_GENES = 6000
MAX_PCT_MT = 20
N_HVG = 2000
N_PCS = 30
RESOLUTION = 0.8
dirs = fair.setup_step("02_qc_cluster")

# %%
# 步骤：读取合并后的原始数据
# 输入：IN_FILE
# 输出：adata
adata = sc.read_h5ad(IN_FILE)
adata

# %%
# 步骤：计算 QC 指标（线粒体、核糖体基因比例）
# 输入：adata
# 输出：adata.obs 新增 n_genes_by_counts、total_counts、pct_counts_mt、pct_counts_ribo
adata.var["mt"] = adata.var_names.str.startswith(MT_PREFIX)
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL", "Rps", "Rpl"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)

# %%
# 步骤：QC 指标小提琴图（按样本），用来确定阈值
# 输入：adata.obs
# 输出：figures/02_qc_violin_by_sample.pdf/.png
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], groupby="sample",
             multi_panel=True, show=False)
fair.save_fig(dirs, "qc_violin_by_sample")

# %%
# 步骤：总 counts vs 基因数散点图（按线粒体比例着色）
# 输入：adata.obs
# 输出：figures/02_qc_scatter_counts_genes.pdf/.png
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt", show=False)
fair.save_fig(dirs, "qc_scatter_counts_genes")

# %%
# 步骤：按阈值过滤细胞和基因
# 输入：adata
# 输出：adata（过滤后）
cells_before = adata.obs["sample"].value_counts()
sc.pp.filter_cells(adata, min_genes=MIN_GENES)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[(adata.obs.n_genes_by_counts < MAX_GENES) & (adata.obs.pct_counts_mt < MAX_PCT_MT)].copy()
adata

# %%
# 步骤：记录每个样本过滤前后的细胞数
# 输入：cells_before、adata
# 输出：tables/02_qc_cells_kept.csv
qc_kept = pd.DataFrame({"before": cells_before, "after": adata.obs["sample"].value_counts()})
qc_kept.to_csv(fair.out(dirs, "tables", "qc_cells_kept.csv"))
qc_kept

# %%
# 步骤：按样本预测双细胞（Scrublet）并去除
# 输入：adata（原始 counts）
# 输出：adata（去掉 predicted_doublet == True 的细胞）；tables/02_doublets_per_sample.csv
sc.pp.scrublet(adata, batch_key="sample", random_state=0)
doublets = pd.crosstab(adata.obs["sample"], adata.obs["predicted_doublet"])
doublets.to_csv(fair.out(dirs, "tables", "doublets_per_sample.csv"))
adata = adata[~adata.obs["predicted_doublet"]].copy()
doublets

# %%
# 步骤：保存原始 counts，再做总量标准化和 log1p
# 输入：adata（原始 counts）
# 输出：adata.X = log1p 标准化；layers['counts']、layers['log1p']
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log1p"] = adata.X.copy()

# %%
# 步骤：用原始 counts 选高变基因（seurat_v3，按样本计算）
# 输入：adata.layers['counts']
# 输出：adata.var['highly_variable']；figures/02_hvg.pdf/.png
sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat_v3", layer="counts", batch_key="sample")
sc.pl.highly_variable_genes(adata, show=False)
fair.save_fig(dirs, "hvg")

# %%
# 步骤：只在高变基因上做 scale 和 PCA（不改动 adata.X）
# 输入：adata
# 输出：adata.obsm['X_pca']
hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(hvg, max_value=10)
sc.pp.pca(hvg, n_comps=50, random_state=0)
adata.obsm["X_pca"] = hvg.obsm["X_pca"]
adata.uns["pca"] = hvg.uns["pca"]

# %%
# 步骤：PCA 方差贡献图，确定 N_PCS
# 输入：adata.uns['pca']
# 输出：figures/02_pca_variance_ratio.pdf/.png
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
fair.save_fig(dirs, "pca_variance_ratio")

# %%
# 步骤：邻居图、Leiden 聚类、UMAP
# 输入：adata.obsm['X_pca']
# 输出：adata.obs['leiden']、adata.obsm['X_umap']；figures/02_umap_leiden_sample.pdf/.png
sc.pp.neighbors(adata, n_pcs=N_PCS, random_state=0)
sc.tl.leiden(adata, resolution=RESOLUTION, flavor="igraph", n_iterations=2, directed=False, random_state=0)
sc.tl.umap(adata, random_state=0)
sc.pl.umap(adata, color=["leiden", "sample", "pct_counts_mt"], show=False)
fair.save_fig(dirs, "umap_leiden_sample")

# %%
# 步骤：每个 leiden 群的 marker 基因（Wilcoxon）
# 输入：adata（X = log1p）
# 输出：tables/02_markers_leiden.csv；figures/02_markers_dotplot_leiden.pdf/.png
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=False)
sc.get.rank_genes_groups_df(adata, group=None).to_csv(fair.out(dirs, "tables", "markers_leiden.csv"), index=False)
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, dendrogram=False, show=False)
fair.save_fig(dirs, "markers_dotplot_leiden")

# %%
# 步骤：保存处理后的数据
# 输入：adata
# 输出：objects/02_processed.h5ad
adata.write(fair.out(dirs, "objects", "processed.h5ad"))

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
