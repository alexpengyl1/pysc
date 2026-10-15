# 16 MAGIC 表达插补
#
# 【这份代码做什么】
#   单细胞数据里很多基因表达是 0（dropout），会掩盖基因之间的关系。
#   MAGIC 借助相似细胞的信息给指定基因做平滑插补，让基因之间的相关关系（如 EMT 中 VIM 升高、CDH1 降低）显现出来。
# 【在流程中的位置】上游 05 → 本步（终点分析）
# 【运行环境】sc-core（课程 ps21）
# 【输入数据】
#   - results/05_annotate/objects/05_annotated.h5ad（X = log1p）
#   - 顶部 GENES：要插补和展示的基因
# 【输出结果】都在 results/16_magic/
#   - tables/16_magic_imputed_<基因>.csv：只含 GENES 的插补表达（细胞 x 基因）
#   - figures/：16_scatter_before_magic、16_scatter_after_magic
# 【怎么解读】
#   - 插补前散点图：大量点贴在坐标轴上（0 值），看不出关系
#   - 插补后：应出现连续的趋势；如果插补后仍无关系，说明这两个基因本身就不相关
# 【注意】
#   - 插补值只用于看基因关系、趋势和展示，不能用于差异分析或统计检验

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/16_magic/）
import scanpy as sc
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
GENES = ["VIM", "CDH1", "ZEB1"]    # 前两个作 x / y，第三个着色
KNN = 5        # 越小越局部、平滑越弱
T = 3          # 扩散步数，越大越平滑
SEED = 0
TAG = "_".join(GENES).lower()
dirs = fair.setup_step("16_magic")

# %%
# 步骤：读取数据
# 输入：IN_FILE
# 输出：adata
adata = sc.read_h5ad(IN_FILE)

# %%
# 步骤：MAGIC 插补指定基因，保存插补值
# 输入：adata.X
# 输出：adata_magic（细胞 x GENES）；tables/16_magic_imputed_<基因>.csv
adata_magic = sc.external.pp.magic(adata, name_list=GENES, knn=KNN, t=T, random_state=SEED, copy=True)
adata_magic.to_df().to_csv(fair.out(dirs, "tables", f"magic_imputed_{TAG}.csv"))

# %%
# 步骤：插补前的基因关系
# 输入：adata
# 输出：figures/16_scatter_before_magic.pdf/.png
sc.pl.scatter(adata, x=GENES[0], y=GENES[1], color=GENES[2], title="Before MAGIC", show=False)
fair.save_fig(dirs, "scatter_before_magic")

# %%
# 步骤：插补后的基因关系
# 输入：adata_magic
# 输出：figures/16_scatter_after_magic.pdf/.png
sc.pl.scatter(adata_magic, x=GENES[0], y=GENES[1], color=GENES[2], title="After MAGIC", show=False)
fair.save_fig(dirs, "scatter_after_magic")

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
