# 19 Hotspot：在一个细胞类型内找局部共变的基因模块
#
# 【这份代码做什么】
#   在一个细胞类型内部，找出在相似细胞之间表达一致变化的基因，并把这些基因聚成模块。
#   每个模块通常对应一种亚型或细胞状态（如活化、应激、增殖）。
#   流程：取细胞类型 → 重算 PCA → 建 KNN 图 → 基因自相关 → 显著基因的局部相关 → 聚成模块 → 模块得分 → UMAP
# 【在流程中的位置】上游 05 → 本步（终点分析；换 CELLTYPE 可对每个细胞类型重复运行）
# 【运行环境】sc-core（课程 ps24）
# 【输入数据】
#   - results/05_annotate/objects/05_annotated.h5ad 中某一细胞类型（layers['counts']、obs['total_counts']；隐空间用重算的 X_pca，也可换 X_scVI）
# 【输出结果】都在 results/19_hotspot/，文件名带细胞类型
#   - tables/：19_autocorrelation_<细胞类型>.csv、19_modules_<细胞类型>.csv、19_module_scores_<细胞类型>.csv
#   - figures/：19_local_correlations_<细胞类型>、19_umap_module_scores_<细胞类型>
# 【怎么解读】
#   - 自相关表：Z 高、FDR 小 = 基因表达在相邻细胞间一致，是区分细胞状态的信息基因
#   - 局部相关热图：对角线上的方块 = 模块；模块 -1 = 未归入任何模块
#   - modules 表：每个基因属于哪个模块，可拿模块基因做富集分析来命名模块
#   - 模块得分 UMAP：模块活跃的细胞区域，对应亚型或活化状态

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/19_hotspot/）
import scanpy as sc
import hotspot
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
CELLTYPE = "T"
N_PCS = 10
N_NEIGHBORS = 30
N_GENES = 500             # 取自相关最强的前 N 个显著基因做模块
MIN_MODULE_GENES = 15
FDR = 0.05
N_JOBS = 4
SEED = 0
TAG = CELLTYPE.lower().replace(" ", "_")
dirs = fair.setup_step("19_hotspot")

# %%
# 步骤：读取数据，取一个细胞类型，重算 PCA
# 输入：IN_FILE
# 输出：adata（子集，obsm['X_pca']、layers['counts_csc']）
adata = sc.read_h5ad(IN_FILE)
adata = adata[adata.obs["cell_type"] == CELLTYPE].copy()
sc.pp.pca(adata, n_comps=N_PCS, mask_var="highly_variable", random_state=SEED)
adata.layers["counts_csc"] = adata.layers["counts"].tocsc()

# %%
# 步骤：建 Hotspot 对象（原始 counts，负二项模型）和 KNN 图
# 输入：adata
# 输出：hs
hs = hotspot.Hotspot(adata, layer_key="counts_csc", model="danb", latent_obsm_key="X_pca",
                     umi_counts_obs_key="total_counts")
hs.create_knn_graph(weighted_graph=False, n_neighbors=N_NEIGHBORS)

# %%
# 步骤：每个基因的自相关
# 输入：hs
# 输出：autocorr；tables/19_autocorrelation_<细胞类型>.csv
autocorr = hs.compute_autocorrelations(jobs=N_JOBS)
autocorr.to_csv(fair.out(dirs, "tables", f"autocorrelation_{TAG}.csv"))
autocorr.head(15)

# %%
# 步骤：显著基因之间的局部相关，聚成模块
# 输入：autocorr
# 输出：modules；tables/19_modules_<细胞类型>.csv
genes = autocorr.loc[autocorr.FDR < FDR].sort_values("Z", ascending=False).head(N_GENES).index
hs.compute_local_correlations(genes, jobs=N_JOBS)
modules = hs.create_modules(min_gene_threshold=MIN_MODULE_GENES, core_only=True, fdr_threshold=FDR)
modules.to_csv(fair.out(dirs, "tables", f"modules_{TAG}.csv"))
modules.value_counts()

# %%
# 步骤：局部相关热图
# 输入：hs
# 输出：figures/19_local_correlations_<细胞类型>.pdf/.png
hs.plot_local_correlations(vmin=-12, vmax=12)
fair.save_fig(dirs, f"local_correlations_{TAG}")

# %%
# 步骤：模块得分
# 输入：hs
# 输出：module_scores；tables/19_module_scores_<细胞类型>.csv
module_scores = hs.calculate_module_scores()
module_scores.to_csv(fair.out(dirs, "tables", f"module_scores_{TAG}.csv"))

# %%
# 步骤：模块得分画在 UMAP 上
# 输入：module_scores
# 输出：figures/19_umap_module_scores_<细胞类型>.pdf/.png
for module in module_scores.columns:
    adata.obs[f"Module {module}"] = module_scores[module]
sc.pl.umap(adata, color=[f"Module {m}" for m in module_scores.columns], vmin=-1, vmax=1, show=False)
fair.save_fig(dirs, f"umap_module_scores_{TAG}")

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
