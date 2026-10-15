# 08a 为 pySCENIC 导出表达矩阵
#
# 【这份代码做什么】
#   pySCENIC 只能在 Python 3.7 环境里运行，读不了新版 h5ad，所以先在 sc-core 里把原始 counts 导出成 csv，
#   同时导出每个细胞的细胞类型，供 08b 和 08c 使用。会去掉表达极低的基因，细胞太多时可以按上限抽样。
# 【在流程中的位置】上游 05 → 本步 → 08b_pyscenic → 08c_scenic_plot
# 【运行环境】sc-core（课程 ps13）
# 【输入数据】
#   - results/05_annotate/objects/05_annotated.h5ad（layers['counts']、obs['cell_type']）
# 【输出结果】都在 results/08_scenic/（08a/08b/08c 共用这个文件夹）
#   - objects/08_expr_counts.csv：细胞 x 基因，原始 counts（08b 的输入）
#   - tables/08_cell_info.csv：细胞 → cell_type（08c 用）
#   - logs/run_info_08a.json、pip_freeze_08a.txt
# 【怎么解读】
#   - 这一步只是格式转换；看打印的矩阵大小，确认细胞数和基因数符合预期
# 【注意】
#   - 抽样上限 MAX_CELLS 只为加快速度；正式分析尽量用全部细胞（课程演示只用了 100 个细胞，不能用于结论）

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/08_scenic/）
import scanpy as sc
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
MIN_CELL_FRACTION = 0.01   # 基因至少在 1% 的细胞里表达
MAX_CELLS = 20000          # GRNBoost2 耗时随细胞数增加；超过就随机抽样
SEED = 0
dirs = fair.setup_step("08_scenic")

# %%
# 步骤：读取数据，去掉表达极低的基因，必要时抽样
# 输入：IN_FILE
# 输出：adata
adata = sc.read_h5ad(IN_FILE)
sc.pp.filter_genes(adata, min_cells=int(MIN_CELL_FRACTION * adata.n_obs))
if adata.n_obs > MAX_CELLS:
    sc.pp.subsample(adata, n_obs=MAX_CELLS, random_state=SEED)
adata

# %%
# 步骤：导出原始 counts 矩阵和细胞注释
# 输入：adata.layers['counts']、adata.obs
# 输出：objects/08_expr_counts.csv、tables/08_cell_info.csv
adata.to_df(layer="counts").to_csv(fair.out(dirs, "objects", "expr_counts.csv"))
adata.obs[["cell_type"]].to_csv(fair.out(dirs, "tables", "cell_info.csv"))

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info_08a.json、logs/pip_freeze_08a.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()}, script="08a")
