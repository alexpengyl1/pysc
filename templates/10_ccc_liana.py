# 10 细胞通讯：LIANA+ 多方法共识排名
#
# 【这份代码做什么】
#   同时运行多种细胞通讯方法（CellPhoneDB、CellChat、NATMI、Connectome、log2FC 等），
#   把各方法的结果合并成一个共识排名，减少单一方法的偏差，再筛选特异且强的配体-受体互作并作图。
# 【在流程中的位置】上游 05 → 本步（终点分析；可与 09 CellPhoneDB 相互印证）
# 【运行环境】liana（课程 ps15）
# 【输入数据】
#   - results/05_annotate/objects/05_annotated.h5ad（X = log1p；obs['cell_type']）
#   - 顶部 SOURCES / TARGETS：要重点看的发送细胞和接收细胞
# 【输出结果】都在 results/10_liana/
#   - tables/：10_liana_rank_aggregate.csv（全部结果）、10_liana_top_specific.csv（specificity_rank <= 0.05）
#   - figures/：10_dotplot_sources_targets、10_circle_plot_counts
# 【怎么解读】
#   - magnitude_rank：互作强度的共识排名，越小越强
#   - specificity_rank：这个互作对这对细胞类型的特异性，越小越特异
#   - 一般先按 specificity_rank <= 0.05 过滤，再按 magnitude_rank 排序
#   - 网络图：线越粗 = 两类细胞之间特异互作越多
# 【注意】
#   - 人用 RESOURCE = "consensus"，小鼠用 "mouseconsensus"

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/10_liana/）
import scanpy as sc
import liana as li
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
CELLTYPE_KEY = "cell_type"
RESOURCE = "consensus"            # 小鼠："mouseconsensus"
EXPR_PROP = 0.1
SPECIFICITY_CUTOFF = 0.05
N_PERMS = 1000
SEED = 0
SOURCES = ["Fibroblast", "Myeloid"]
TARGETS = ["T", "Epithelial"]
dirs = fair.setup_step("10_liana")

# %%
# 步骤：读取数据
# 输入：IN_FILE
# 输出：adata
adata = sc.read_h5ad(IN_FILE)
adata

# %%
# 步骤：多方法共识排名（X 是 log1p，所以 use_raw=False）
# 输入：adata
# 输出：adata.uns['liana_res']；tables/10_liana_rank_aggregate.csv
li.mt.rank_aggregate(adata, groupby=CELLTYPE_KEY, resource_name=RESOURCE, expr_prop=EXPR_PROP,
                     use_raw=False, n_perms=N_PERMS, seed=SEED, verbose=True)
liana_res = adata.uns["liana_res"]
liana_res.to_csv(fair.out(dirs, "tables", "liana_rank_aggregate.csv"), index=False)
liana_res.head()

# %%
# 步骤：只保留特异的互作，按强度排序
# 输入：liana_res
# 输出：tables/10_liana_top_specific.csv
top = liana_res[liana_res["specificity_rank"] <= SPECIFICITY_CUTOFF].sort_values("magnitude_rank")
top.to_csv(fair.out(dirs, "tables", "liana_top_specific.csv"), index=False)
top.head()

# %%
# 步骤：指定来源 / 接收细胞的互作点图
# 输入：adata.uns['liana_res']
# 输出：figures/10_dotplot_sources_targets.pdf/.png
dotplot = li.pl.dotplot(adata=adata, colour="magnitude_rank", size="specificity_rank",
                        inverse_colour=True, inverse_size=True,
                        source_labels=SOURCES, target_labels=TARGETS,
                        filter_fun=lambda x: x["specificity_rank"] <= SPECIFICITY_CUTOFF,
                        top_n=20, orderby="magnitude_rank", orderby_ascending=True, figure_size=(10, 7))
fair.save_fig(dirs, "dotplot_sources_targets", fig=dotplot)

# %%
# 步骤：细胞类型之间的特异互作数量网络图
# 输入：adata.uns['liana_res']
# 输出：figures/10_circle_plot_counts.pdf/.png
li.pl.circle_plot(adata, groupby=CELLTYPE_KEY, score_key="magnitude_rank", inverse_score=True,
                  source_labels=SOURCES, filter_fun=lambda x: x["specificity_rank"] <= SPECIFICITY_CUTOFF,
                  pivot_mode="counts", figure_size=(8, 8))
fair.save_fig(dirs, "circle_plot_counts")

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
