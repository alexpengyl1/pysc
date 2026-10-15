# 24 CellOracle：按细胞群构建基因调控网络（GRN）并比较网络中心性
#
# 【这份代码做什么】
#   结合先验的 TF 结合位点网络（base GRN）和单细胞表达数据，为每个细胞群分别构建 TF → 靶基因 的调控网络，
#   再计算每个基因在网络中的重要性（中心性），比较不同细胞群之间关键调控因子的变化。
#   流程：24a 准备输入（sc-core）→ 在 Docker 中：建 Oracle → PCA 和 KNN 插补 → 每个群拟合 GRN → 过滤边、算网络分数 → 排名和比较图
# 【在流程中的位置】上游 05 → 本步（终点分析）
# 【运行环境】24a 块在 sc-core 里运行；其余在 Docker 容器 kenjikamimoto126/celloracle_ubuntu:0.18.0 的 Jupyter 里运行
#   （启动方法见 SKILL.md；容器挂载整个项目文件夹，fair_utils.py 也在里面）（课程 ps30）
# 【输入数据】
#   - results/05_annotate/objects/05_annotated.h5ad（layers['counts']、obs['cell_type']、obsm['X_umap']）
#   - 先验 base GRN：CellOracle 自带（人或小鼠）
# 【输出结果】都在 results/24_celloracle/
#   - objects/：24_input.h5ad（24a 生成）、24_oracle.celloracle.oracle、24_links.celloracle.links
#   - tables/：24_grn_<细胞群>.csv（过滤后的 TF → 靶基因 边表）、24_network_scores.csv
#   - figures/：24_pca_cumulative_variance、24_score_comparison_<群1>_vs_<群2>、24_network_score_<TF>，以及 CellOracle 自带函数存的排名图
# 【怎么解读】
#   - degree_centrality_all：连接数多 = hub 基因
#   - eigenvector_centrality：连接到其他重要节点的基因
#   - betweenness_centrality：处在信息流“桥梁”位置的基因
#   - score_comparison 图：偏离对角线的 TF = 在两个群中重要性差别大，是候选的状态特异调控因子
#   - network_score_<TF>：某个 TF 的网络重要性在各群间的变化
# 【注意】
#   - 课程只讲了网络构建；扰动模拟（simulate_shift）见 CellOracle 官方教程
#   - 每个群拟合 GRN 约需半小时

# %%
# 24a（sc-core 环境）
# 步骤：准备 CellOracle 输入：≤3 万细胞、3000 个高变基因、原始 counts 放在 layers['raw_count']
# 输入：results/05_annotate/objects/05_annotated.h5ad
# 输出：results/24_celloracle/objects/24_input.h5ad
# import scanpy as sc
# import fair_utils as fair
# dirs = fair.setup_step("24_celloracle")
# adata = sc.read_h5ad("results/05_annotate/objects/05_annotated.h5ad")
# if adata.n_obs > 30000:
#     sc.pp.subsample(adata, n_obs=30000, random_state=0)
# sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3", layer="counts")
# adata = adata[:, adata.var.highly_variable].copy()
# adata.layers["raw_count"] = adata.layers["counts"].copy()
# adata.write(fair.out(dirs, "objects", "input.h5ad"))

# %%
# 以下在 Docker 容器的 Jupyter 中运行（工作目录 = 项目文件夹）
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/24_celloracle/）
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import celloracle as co
import fair_utils as fair

IN_FILE = "results/24_celloracle/objects/24_input.h5ad"
CLUSTER_KEY = "cell_type"
EMBEDDING = "X_umap"
SPECIES = "mouse"                 # mouse：小鼠 scATAC 图谱 base GRN；human：人启动子 base GRN
ALPHA = 10                        # GRN 正则化强度，越大越稀疏
P_CUTOFF = 0.001
MAX_EDGES = 2000                  # 每个群最多保留的边数
GOI = "Gata2"                     # 想看的转录因子
COMPARE = ("Stem", "Enterocyte")  # 比较中心性的两个群
dirs = fair.setup_step("24_celloracle")

# %%
# 步骤：读取数据和先验 base GRN
# 输入：IN_FILE
# 输出：adata、base_grn
adata = sc.read_h5ad(IN_FILE)
if SPECIES == "mouse":
    base_grn = co.data.load_mouse_scATAC_atlas_base_GRN()
else:
    base_grn = co.data.load_human_promoter_base_GRN()
base_grn.head()

# %%
# 步骤：建 Oracle 对象，导入原始 counts 和先验网络
# 输入：adata、base_grn
# 输出：oracle
adata.X = adata.layers["raw_count"].copy()
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name=CLUSTER_KEY, embedding_name=EMBEDDING)
oracle.import_TF_data(TF_info_matrix=base_grn)

# %%
# 步骤：PCA，按累计方差拐点选主成分数
# 输入：oracle
# 输出：n_comps；figures/24_pca_cumulative_variance.pdf/.png
oracle.perform_PCA()
cumulative = np.cumsum(oracle.pca.explained_variance_ratio_)[:100]
n_comps = min(np.where(np.diff(np.diff(cumulative) > 0.002))[0][0], 50)
plt.plot(cumulative)
plt.axvline(n_comps, color="black")
fair.save_fig(dirs, "pca_cumulative_variance")
n_comps

# %%
# 步骤：KNN 插补（k = 细胞数的 2.5%），保存 oracle
# 输入：oracle、n_comps
# 输出：objects/24_oracle.celloracle.oracle
k = int(0.025 * oracle.adata.n_obs)
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k * 8, b_maxl=k * 4, n_jobs=4)
oracle.to_hdf5(str(fair.out(dirs, "objects", "oracle.celloracle.oracle")))

# %%
# 步骤：每个细胞群拟合 GRN（耗时，约半小时）
# 输入：oracle
# 输出：links（每个群一张边表）
links = oracle.get_links(cluster_name_for_GRN_unit=CLUSTER_KEY, alpha=ALPHA, verbose_level=10)

# %%
# 步骤：过滤边，计算网络分数，保存 links
# 输入：links
# 输出：objects/24_links.celloracle.links
links.filter_links(p=P_CUTOFF, weight="coef_abs", threshold_number=MAX_EDGES)
links.get_network_score()
links.to_hdf5(file_path=str(fair.out(dirs, "objects", "links.celloracle.links")))

# %%
# 步骤：保存每个群的边表和全部网络分数
# 输入：links
# 输出：tables/24_grn_<细胞群>.csv、24_network_scores.csv
for cluster, edges in links.filtered_links.items():
    edges.to_csv(fair.out(dirs, "tables", f"grn_{cluster}.csv"), index=False)
links.merged_score.to_csv(fair.out(dirs, "tables", "network_scores.csv"))

# %%
# 步骤：各群中心性最高的基因排名图（CellOracle 自带函数，直接存到 figures/）
# 输入：links
# 输出：figures/ 下的排名图
for cluster in links.links_dict.keys():
    links.plot_scores_as_rank(cluster=cluster, n_gene=30, save=str(dirs["figures"]))

# %%
# 步骤：两个群之间比较特征向量中心性
# 输入：links、COMPARE
# 输出：figures/24_score_comparison_<群1>_vs_<群2>.pdf/.png
links.plot_score_comparison_2D(value="eigenvector_centrality", cluster1=COMPARE[0], cluster2=COMPARE[1], percentile=98)
fair.save_fig(dirs, f"score_comparison_{COMPARE[0]}_vs_{COMPARE[1]}")

# %%
# 步骤：某个 TF 的网络分数在各群间的变化
# 输入：links、GOI
# 输出：figures/24_network_score_<TF>.pdf/.png
links.plot_score_per_cluster(goi=GOI)
fair.save_fig(dirs, f"network_score_{GOI}")

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
