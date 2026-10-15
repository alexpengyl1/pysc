# 21a 为 scFEA 导出表达矩阵
#
# 【这份代码做什么】
#   scFEA 是独立的命令行工具（单独的 scFEA 环境），需要 基因 x 细胞 的 csv 表达矩阵。
#   这一步在 sc-core 里把 log 标准化表达导出成这个格式，同时导出每个细胞的细胞类型，细胞太多时可以按上限抽样。
#   （课程 ps27 是从 Seurat 的 data slot 导出，这里用 scanpy 做等价导出）
# 【在流程中的位置】上游 05 → 本步 → 21_metabolism_scfea.sh
# 【运行环境】sc-core
# 【输入数据】
#   - results/05_annotate/objects/05_annotated.h5ad（layers['log1p']、obs['cell_type']）
# 【输出结果】都在 results/21_scfea/（21a 和 21 共用）
#   - objects/21_expr.csv：基因 x 细胞（21_metabolism_scfea.sh 的输入）
#   - tables/21_cell_info.csv：细胞 → cell_type（之后按细胞类型汇总通量时用）
#   - logs/run_info_21a.json、pip_freeze_21a.txt
# 【怎么解读】
#   - 这一步只是格式转换；看打印的矩阵大小，确认基因数和细胞数符合预期
# 【注意】
#   - scFEA 每次都要训练模型，细胞越多越慢；抽样数 N_CELLS 按电脑性能调整

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/21_scfea/）
import scanpy as sc
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
N_CELLS = 2000     # scFEA 每次运行都训练模型，先抽样
SEED = 0
dirs = fair.setup_step("21_scfea")

# %%
# 步骤：读取数据并抽样
# 输入：IN_FILE
# 输出：adata（N_CELLS 个细胞）
adata = sc.read_h5ad(IN_FILE)
sc.pp.subsample(adata, n_obs=min(N_CELLS, adata.n_obs), random_state=SEED)

# %%
# 步骤：导出 基因 x 细胞 的 log 标准化矩阵和细胞注释
# 输入：adata
# 输出：objects/21_expr.csv、tables/21_cell_info.csv
adata.to_df(layer="log1p").T.to_csv(fair.out(dirs, "objects", "expr.csv"))
adata.obs[["cell_type"]].to_csv(fair.out(dirs, "tables", "cell_info.csv"))

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info_21a.json、logs/pip_freeze_21a.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()}, script="21a")
