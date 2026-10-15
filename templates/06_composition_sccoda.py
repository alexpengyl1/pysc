# 06 细胞组成变化分析（scCODA）
#
# 【这份代码做什么】
#   检验不同条件（如照射 vs 对照）之间各细胞类型的比例是否真的改变。
#   细胞比例是成分数据（一种细胞变多，其他的比例自然变少），scCODA 用贝叶斯模型并指定一个参考细胞类型来处理这个问题。
#   流程：汇总 样本 x 细胞类型 计数 → 堆叠柱状图和箱线图 → 选参考细胞类型 → 建模采样 → 输出可信变化
# 【在流程中的位置】上游 05 → 本步（终点分析）
# 【运行环境】sc-sccoda（TensorFlow；课程 ps11）
# 【输入数据】
#   - results/05_annotate/objects/05_annotated.h5ad（obs['cell_type']、obs['sample']）
#   - data/sample_info.csv（sample, condition；每个样本一行）
# 【输出结果】都在 results/06_sccoda/
#   - tables/：06_counts_sample_celltype.csv（计数表）、06_effects.csv、06_intercepts.csv、06_credible_effects.csv
#   - figures/：06_stacked_barplot_condition、06_boxplot_condition、06_reference_selection
#   - objects/06_sccoda_result.pkl：模型结果，可重新加载
# 【怎么解读】
#   - 堆叠柱状图 / 箱线图：直观看各组比例差异，箱线图里每个点是一个样本
#   - effects 表：Final Parameter 不为 0 = 可信变化；log2-fold change 是相对对照组的倍数
#   - credible_effects：最终结论以它为准；其他细胞类型的“下降”可能只是比例被挤占
#   - reference_selection 图：参考细胞类型要选在所有样本里都存在、比例稳定的（出现率高、离散度低），结论依赖这个选择
# 【注意】
#   - 每个条件至少 2 个生物学重复样本；统计单位是样本，不是细胞
#   - 默认 FDR 0.05 很严格；改成 0.2 更灵敏，但假阳性更多

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/06_sccoda/）
import pandas as pd
import scanpy as sc
from sccoda.util import cell_composition_data as dat
from sccoda.util import comp_ana as mod
from sccoda.util import data_visualization as viz
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
SAMPLE_INFO = "data/sample_info.csv"
CELLTYPE_KEY = "cell_type"
SAMPLE_KEY = "sample"
CONDITION_KEY = "condition"
REFERENCE = "automatic"          # 或指定细胞类型名
FDR = 0.05
dirs = fair.setup_step("06_sccoda")

# %%
# 步骤：读取注释好的数据和样本分组表
# 输入：IN_FILE、SAMPLE_INFO
# 输出：adata、sample_info（index = 样本名）
adata = sc.read_h5ad(IN_FILE)
sample_info = pd.read_csv(SAMPLE_INFO, index_col="sample")
sample_info

# %%
# 步骤：汇总成 样本 x 细胞类型 计数表并保存
# 输入：adata.obs、sample_info
# 输出：counts（AnnData：行 = 样本，列 = 细胞类型）；tables/06_counts_sample_celltype.csv
counts = dat.from_scanpy(adata, cell_type_identifier=CELLTYPE_KEY, sample_identifier=SAMPLE_KEY,
                         covariate_df=sample_info)
counts.to_df().to_csv(fair.out(dirs, "tables", "counts_sample_celltype.csv"))
counts

# %%
# 步骤：各条件的细胞组成堆叠柱状图
# 输入：counts
# 输出：figures/06_stacked_barplot_condition.pdf/.png
viz.stacked_barplot(counts, feature_name=CONDITION_KEY)
fair.save_fig(dirs, "stacked_barplot_condition")

# %%
# 步骤：各细胞类型比例按条件的箱线图
# 输入：counts
# 输出：figures/06_boxplot_condition.pdf/.png
viz.boxplots(counts, feature_name=CONDITION_KEY, add_dots=True)
fair.save_fig(dirs, "boxplot_condition")

# %%
# 步骤：参考细胞类型选择图（出现率高、离散度低的适合做参考）
# 输入：counts
# 输出：figures/06_reference_selection.pdf/.png
viz.rel_abundance_dispersion_plot(data=counts, abundant_threshold=0.9)
fair.save_fig(dirs, "reference_selection")

# %%
# 步骤：建模并 HMC 采样
# 输入：counts
# 输出：result（CAResult）
model = mod.CompositionalAnalysis(counts, formula=CONDITION_KEY, reference_cell_type=REFERENCE)
result = model.sample_hmc()
result.set_fdr(est_fdr=FDR)
result.summary()

# %%
# 步骤：保存效应、截距和可信效应表
# 输入：result
# 输出：tables/06_effects.csv、06_intercepts.csv、06_credible_effects.csv
result.effect_df.to_csv(fair.out(dirs, "tables", "effects.csv"))
result.intercept_df.to_csv(fair.out(dirs, "tables", "intercepts.csv"))
result.credible_effects().to_csv(fair.out(dirs, "tables", "credible_effects.csv"))

# %%
# 步骤：保存结果对象
# 输入：result
# 输出：objects/06_sccoda_result.pkl
result.save(str(fair.out(dirs, "objects", "sccoda_result.pkl")))

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
