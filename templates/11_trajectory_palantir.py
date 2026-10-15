# 11 Palantir：拟时序、分化潜能和终末命运概率
#
# 【这份代码做什么】
#   对一个分化谱系（例如干细胞 → 各种成熟细胞）推断每个细胞的分化进程（拟时序）、
#   还有多大可能分化成不同方向（熵），以及走向每个终末状态的概率，并画出基因沿各分支的变化趋势。
#   流程：取谱系子集 → diffusion map → MAGIC 插补 → 确认起点和终点 → 运行 Palantir → 分支轨迹 → 基因趋势
# 【在流程中的位置】上游 05 → 本步 → 下游 14_fate_cellrank（可选）
# 【运行环境】sc-palantir（课程 ps16）
# 【输入数据】
#   - results/05_annotate/objects/05_annotated.h5ad 中要分析的谱系子集（obs['cell_type']、obsm['X_umap']）
#   - 顶部 START_CELL：起始细胞 barcode（例如干/祖细胞群里干性 marker 最高的细胞）
#   - 顶部 TERMINAL_STATES：各终末状态的代表细胞 barcode
# 【输出结果】都在 results/11_palantir/
#   - objects/11_palantir.h5ad：obs 有 palantir_pseudotime、palantir_entropy；obsm 有 palantir_fate_probabilities；layers['MAGIC_imputed_data']（14 会用到）
#   - tables/：11_pseudotime_entropy.csv、11_fate_probabilities.csv
#   - figures/：11_start_terminal_cells、11_palantir_results、11_branch_selection、11_trajectories、11_gene_trends、11_gene_trend_heatmaps
# 【怎么解读】
#   - start_terminal_cells 图：先确认起点和终点落在正确的细胞群上，选错了后面全错
#   - pseudotime：从起始细胞出发的分化进程，0 → 1
#   - entropy：高 = 多潜能、命运未定；沿分化方向应逐渐下降
#   - fate probabilities：每个细胞走向各终末状态的概率，和为 1
#   - gene trends：基因沿各分支随拟时序的变化，找分支特异的上调/下调基因
# 【注意】
#   - 只放同一谱系的细胞，不相关的细胞类型会干扰轨迹
#   - MAGIC 插补值只用于画趋势，不用于统计检验

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/11_palantir/）
import pandas as pd
import scanpy as sc
import palantir
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
LINEAGE_TYPES = ["Stem", "TA", "Enterocyte", "Goblet"]   # 要分析的谱系细胞类型
START_CELL = "AAACCCAAGXXXXXXX-1-S1"                     # 起始细胞 barcode
TERMINAL_STATES = {                                       # 代表细胞 barcode: 终末状态名
    "AAACGAAXXXXXXX-1-S1": "Enterocyte",
    "AAAGTCCXXXXXXX-1-S2": "Goblet",
}
N_WAYPOINTS = 500
GENES = ["LGR5", "MKI67", "FABP2", "MUC2"]
SEED = 0
dirs = fair.setup_step("11_palantir")

# %%
# 步骤：读取数据，取出谱系子集，重算 PCA
# 输入：IN_FILE
# 输出：adata（谱系子集）
adata = sc.read_h5ad(IN_FILE)
adata = adata[adata.obs["cell_type"].isin(LINEAGE_TYPES)].copy()
sc.pp.pca(adata, mask_var="highly_variable", random_state=SEED)
adata

# %%
# 步骤：diffusion map 和多尺度空间
# 输入：adata.obsm['X_pca']
# 输出：adata.obsm['DM_EigenVectors_multiscaled'] 等
palantir.utils.run_diffusion_maps(adata, n_components=5)
palantir.utils.determine_multiscale_space(adata)

# %%
# 步骤：MAGIC 插补（只用于基因趋势和展示）
# 输入：adata
# 输出：adata.layers['MAGIC_imputed_data']
palantir.utils.run_magic_imputation(adata)

# %%
# 步骤：在 UMAP 上确认起始细胞和终末细胞的位置
# 输入：START_CELL、TERMINAL_STATES
# 输出：figures/11_start_terminal_cells.pdf/.png
terminal_states = pd.Series(TERMINAL_STATES)
palantir.plot.highlight_cells_on_umap(adata, pd.concat([terminal_states, pd.Series({START_CELL: "start"})]))
fair.save_fig(dirs, "start_terminal_cells")

# %%
# 步骤：运行 Palantir
# 输入：adata、START_CELL、terminal_states
# 输出：adata.obs['palantir_pseudotime']、['palantir_entropy']、obsm['palantir_fate_probabilities']
palantir.core.run_palantir(adata, START_CELL, num_waypoints=N_WAYPOINTS, terminal_states=terminal_states, seed=SEED)

# %%
# 步骤：拟时序、熵、命运概率总览图
# 输入：adata
# 输出：figures/11_palantir_results.pdf/.png
palantir.plot.plot_palantir_results(adata, s=3)
fair.save_fig(dirs, "palantir_results")

# %%
# 步骤：保存每个细胞的拟时序、熵和命运概率
# 输入：adata
# 输出：tables/11_pseudotime_entropy.csv、11_fate_probabilities.csv
adata.obs[["cell_type", "palantir_pseudotime", "palantir_entropy"]].to_csv(fair.out(dirs, "tables", "pseudotime_entropy.csv"))
pd.DataFrame(adata.obsm["palantir_fate_probabilities"], index=adata.obs_names).to_csv(fair.out(dirs, "tables", "fate_probabilities.csv"))

# %%
# 步骤：选出每个分支的代表细胞
# 输入：adata
# 输出：adata.obsm['branch_masks']；figures/11_branch_selection.pdf/.png
palantir.presults.select_branch_cells(adata, q=0.01, eps=0.01)
palantir.plot.plot_branch_selection(adata)
fair.save_fig(dirs, "branch_selection")

# %%
# 步骤：各分支轨迹
# 输入：adata
# 输出：figures/11_trajectories.pdf/.png
palantir.plot.plot_trajectories(adata, pseudotime_interval=(0, 0.9))
fair.save_fig(dirs, "trajectories")

# %%
# 步骤：基因沿各分支的表达趋势（曲线图）
# 输入：adata.layers['MAGIC_imputed_data']
# 输出：figures/11_gene_trends.pdf/.png
palantir.presults.compute_gene_trends(adata, expression_key="MAGIC_imputed_data")
palantir.plot.plot_gene_trends(adata, GENES)
fair.save_fig(dirs, "gene_trends")

# %%
# 步骤：基因趋势热图
# 输入：adata
# 输出：figures/11_gene_trend_heatmaps.pdf/.png
palantir.plot.plot_gene_trend_heatmaps(adata, GENES)
fair.save_fig(dirs, "gene_trend_heatmaps")

# %%
# 步骤：保存
# 输入：adata
# 输出：objects/11_palantir.h5ad
adata.write(fair.out(dirs, "objects", "palantir.h5ad"))

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
