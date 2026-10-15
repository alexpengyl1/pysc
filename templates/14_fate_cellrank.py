# 14 CellRank：终末状态、命运概率和谱系驱动基因
#
# 【这份代码做什么】
#   在 Palantir 拟时序（11）或 RNA velocity（13）的基础上构建细胞间的转移矩阵，
#   自动找出稳定的宏观状态和终末状态，计算每个细胞走向各终末状态的概率，并找出驱动每条谱系的基因。
#   流程：建转移矩阵（方向性核 + 相似性核）→ 随机游走检查方向 → GPCCA 宏观状态 → 终末状态 → 命运概率 → 驱动基因 → 基因趋势
# 【在流程中的位置】上游 11（Palantir）或 13（scVelo）→ 本步（终点分析）
# 【运行环境】用 Palantir 拟时序 → sc-palantir；用 RNA velocity → sc-velo（课程 ps19）
# 【输入数据】
#   - results/11_palantir/objects/11_palantir.h5ad（obs['palantir_pseudotime']）
#   - 或 results/13_scvelo/objects/13_velocity.h5ad（velocity 层；改用代码里注释掉的 VelocityKernel 块）
# 【输出结果】都在 results/14_cellrank/
#   - objects/14_cellrank.h5ad
#   - tables/：14_terminal_states.csv、14_fate_probabilities.csv、14_drivers_<谱系>.csv
#   - figures/：14_random_walks、14_macrostates_all、14_macrostates_terminal、14_fate_probabilities、14_circular_projection、14_gene_trends（CellRank 自带作图函数只存 PDF）
# 【怎么解读】
#   - random walks：从起始群出发的模拟路径应流向预期的终末群，否则先检查输入的拟时序或速度
#   - macrostates：粗粒化出的稳定状态；N_STATES 太小会把谱系混在一起，太大会拆出噪声态
#   - fate probabilities：每个细胞到各终末状态的概率
#   - circular projection：各细胞类型偏向哪个命运
#   - lineage drivers：与某谱系命运概率最相关的基因（corr 高、qval 小），可作为候选调控基因

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/14_cellrank/）
import numpy as np
import pandas as pd
import scanpy as sc
import cellrank as cr
import fair_utils as fair

IN_FILE = "results/11_palantir/objects/11_palantir.h5ad"
TIME_KEY = "palantir_pseudotime"
CLUSTER_KEY = "cell_type"
START_CLUSTER = "Stem"
KERNEL_WEIGHT = 0.8        # 方向性核（拟时序或 velocity）的权重，其余给表达相似性核
N_STATES = 6
GENES = ["LGR5", "MKI67", "FABP2", "MUC2"]
SEED = 0
dirs = fair.setup_step("14_cellrank")
cr.settings.figdir = str(dirs["figures"])

# %%
# 步骤：读取数据，重算邻居图
# 输入：IN_FILE
# 输出：adata
adata = sc.read_h5ad(IN_FILE)
sc.pp.neighbors(adata, random_state=SEED)
adata

# %%
# 步骤：构建转移矩阵（方向性核 + 表达相似性核）
# 输入：adata.obs[TIME_KEY]、adata.obsp['connectivities']
# 输出：pk、kernel
pk = cr.kernels.PseudotimeKernel(adata, time_key=TIME_KEY).compute_transition_matrix()
ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
kernel = KERNEL_WEIGHT * pk + (1 - KERNEL_WEIGHT) * ck
# 用 RNA velocity 时改为：
# vk = cr.kernels.VelocityKernel(adata).compute_transition_matrix()
# kernel = KERNEL_WEIGHT * vk + (1 - KERNEL_WEIGHT) * ck

# %%
# 步骤：模拟随机游走，检查方向是否合理
# 输入：pk
# 输出：figures/14_random_walks.pdf
pk.plot_random_walks(n_sims=100, start_ixs={CLUSTER_KEY: START_CLUSTER}, basis="umap", seed=SEED,
                     save="14_random_walks.pdf")

# %%
# 步骤：GPCCA 粗粒化，画全部宏观状态
# 输入：kernel
# 输出：g（估计器）；figures/14_macrostates_all.pdf
g = cr.estimators.GPCCA(kernel)
g.fit(n_states=N_STATES, cluster_key=CLUSTER_KEY)
g.plot_macrostates(which="all", save="14_macrostates_all.pdf")

# %%
# 步骤：预测终末状态，保存每个细胞的终末状态标签
# 输入：g
# 输出：figures/14_macrostates_terminal.pdf；tables/14_terminal_states.csv
g.predict_terminal_states()
g.plot_macrostates(which="terminal", save="14_macrostates_terminal.pdf")
g.terminal_states.to_csv(fair.out(dirs, "tables", "terminal_states.csv"))

# %%
# 步骤：命运概率（图 + 表）
# 输入：g
# 输出：figures/14_fate_probabilities.pdf；tables/14_fate_probabilities.csv
g.compute_fate_probabilities()
g.plot_fate_probabilities(same_plot=False, save="14_fate_probabilities.pdf")
fate = pd.DataFrame(np.asarray(g.fate_probabilities), index=adata.obs_names, columns=g.fate_probabilities.names)
fate.to_csv(fair.out(dirs, "tables", "fate_probabilities.csv"))

# %%
# 步骤：各细胞类型的命运偏好（圆形投影）
# 输入：adata（命运概率）
# 输出：figures/14_circular_projection.pdf
cr.pl.circular_projection(adata, keys=CLUSTER_KEY, legend_loc="right", save="14_circular_projection.pdf")

# %%
# 步骤：每条谱系的驱动基因
# 输入：g
# 输出：每条谱系一个 tables/14_drivers_<谱系>.csv
for lineage in g.terminal_states.cat.categories:
    drivers = g.compute_lineage_drivers(lineages=[lineage])
    drivers.to_csv(fair.out(dirs, "tables", f"drivers_{lineage}.csv"))

# %%
# 步骤：基因趋势（pyGAM，不需要 R；要用 R mgcv 时换成 cr.models.GAMR）
# 输入：adata、命运概率
# 输出：figures/14_gene_trends.pdf
model = cr.models.GAM(adata)
cr.pl.gene_trends(adata, model=model, genes=GENES, time_key=TIME_KEY, same_plot=True, hide_cells=True,
                  save="14_gene_trends.pdf")

# %%
# 步骤：保存
# 输入：adata
# 输出：objects/14_cellrank.h5ad
adata.write(fair.out(dirs, "objects", "cellrank.h5ad"))

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
