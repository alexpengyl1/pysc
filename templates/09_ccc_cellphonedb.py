# 09 细胞通讯：CellPhoneDB v5（统计方法）+ ktplotspy 作图
#
# 【这份代码做什么】
#   根据配体和受体在不同细胞类型中的表达，找出细胞类型之间可能存在的信号交流（配体-受体对），
#   用置换检验判断哪些互作是特异的，并统计每对细胞类型之间显著互作的数量。
#   流程：下载数据库 → （小鼠）转成人类同源基因 → 写出输入文件 → 置换检验 → 统计显著互作 → 热图和点图
# 【在流程中的位置】上游 05 → 本步（终点分析；可与 10 LIANA 相互印证）
# 【运行环境】cpdb（课程 ps14）
# 【输入数据】
#   - results/05_annotate/objects/05_annotated.h5ad（X = log1p，不要 scale；obs['cell_type']）
#   - 可选：data/microenvironment.tsv（cell_type, microenvironment），只比较同一微环境内的细胞；data/active_TFs.tsv（cell_type, TF）
# 【输出结果】都在 results/09_cellphonedb/
#   - objects/：09_counts_log1p.h5ad、09_metadata.tsv（CellPhoneDB 的输入文件）
#   - tables/：CellPhoneDB 原始输出 statistical_analysis_pvalues / means / significant_means / interaction_scores / deconvoluted（文件名带 09）、09_n_significant_by_pair.csv
#   - figures/：09_heatmap_n_significant、09_dotplot_<来源>_to_<接收>
# 【怎么解读】
#   - 热图：每对细胞类型之间显著互作的数量，数量多 = 交流频繁（只说明“可能”，不是实验证据）
#   - significant_means：行 = 配体|受体，列 = 细胞A|细胞B；有数值 = 显著，空 = 不显著
#   - interaction_scores：0–100，越高越特异
#   - 点图：点大小 = 显著性，颜色 = 平均表达；挑出与研究问题相关的通路再验证
# 【注意】
#   - 数据库只有人类基因。小鼠数据先用同源表（Ensembl BioMart 导出 mouse_symbol, human_symbol）换成人类基因名，不能简单转大写

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/09_cellphonedb/）
import pandas as pd
import anndata as ad
import ktplotspy as kpy
from cellphonedb.utils import db_utils
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
DB_DIR = "data/cpdb_db"
DB_VERSION = "v5.0.0"
CELLTYPE_KEY = "cell_type"
SPECIES = "human"              # 小鼠改 "mouse"，并准备同源表 ORTHOLOG_FILE
ORTHOLOG_FILE = "data/mouse_human_orthologs.csv"
ITERATIONS = 1000
EXPR_THRESHOLD = 0.1           # 细胞类型内表达比例 < 10% 的基因视为不表达
PVALUE = 0.05
SEED = 42
SOURCE = "Fibroblast"                       # 作图用：配体来源细胞
TARGETS = "Myeloid|T"                       # 作图用：受体细胞，用 | 分隔
GENES = ["TGFB1", "CXCL12", "CSF1"]
dirs = fair.setup_step("09_cellphonedb")

# %%
# 步骤：下载 CellPhoneDB 数据库（只需一次）
# 输入：无
# 输出：data/cpdb_db/cellphonedb.zip
db_utils.download_database(DB_DIR, DB_VERSION)

# %%
# 步骤：读取注释数据
# 输入：IN_FILE
# 输出：adata
adata = ad.read_h5ad(IN_FILE)
adata

# %%
# 步骤（仅小鼠）：基因名换成人类同源基因，去掉没有同源或重复的基因
# 输入：adata、ORTHOLOG_FILE
# 输出：adata（var_names = 人类基因名）
if SPECIES == "mouse":
    orthologs = pd.read_csv(ORTHOLOG_FILE).drop_duplicates("mouse_symbol").set_index("mouse_symbol")
    adata = adata[:, adata.var_names.isin(orthologs.index)].copy()
    adata.var_names = orthologs.loc[adata.var_names, "human_symbol"].values
    adata = adata[:, ~adata.var_names.duplicated()].copy()

# %%
# 步骤：写出 CellPhoneDB 需要的表达矩阵和细胞注释
# 输入：adata
# 输出：objects/09_counts_log1p.h5ad、objects/09_metadata.tsv
counts_file = fair.out(dirs, "objects", "counts_log1p.h5ad")
meta_file = fair.out(dirs, "objects", "metadata.tsv")
adata.write(counts_file)
pd.DataFrame({"barcode_sample": adata.obs_names, "cell_type": adata.obs[CELLTYPE_KEY].values}).to_csv(
    meta_file, sep="\t", index=False)

# %%
# 步骤：运行统计分析（置换检验）
# 输入：数据库 zip、metadata、counts_log1p
# 输出：cpdb_results（字典）；CellPhoneDB 结果表写入 tables/
cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path=f"{DB_DIR}/cellphonedb.zip",
    meta_file_path=str(meta_file),
    counts_file_path=str(counts_file),
    counts_data="hgnc_symbol",
    # microenvs_file_path="data/microenvironment.tsv",   # 可选：只在同一微环境内配对
    # active_tfs_file_path="data/active_TFs.tsv",        # 可选：CellSign
    score_interactions=True,
    iterations=ITERATIONS,
    threshold=EXPR_THRESHOLD,
    threads=4,
    debug_seed=SEED,
    pvalue=PVALUE,
    output_path=str(dirs["tables"]),
    output_suffix="09",
)
list(cpdb_results.keys())

# %%
# 步骤：每对细胞类型的显著互作数量
# 输入：cpdb_results['pvalues']
# 输出：tables/09_n_significant_by_pair.csv
pvals = cpdb_results["pvalues"].set_index("interacting_pair").filter(like="|")
n_sig = (pvals < PVALUE).sum().rename("n_significant").sort_values(ascending=False)
n_sig.to_csv(fair.out(dirs, "tables", "n_significant_by_pair.csv"))
n_sig.head()

# %%
# 步骤：各细胞类型之间显著互作数量热图
# 输入：cpdb_results['pvalues']
# 输出：figures/09_heatmap_n_significant.pdf/.png
grid = kpy.plot_cpdb_heatmap(pvals=cpdb_results["pvalues"], figsize=(6, 6), title="Number of significant interactions")
fair.save_fig(dirs, "heatmap_n_significant", fig=grid)

# %%
# 步骤：指定细胞类型和基因的互作点图
# 输入：adata、cpdb_results
# 输出：figures/09_dotplot_<来源>_to_<接收>.pdf/.png
dotplot = kpy.plot_cpdb(
    adata=adata,
    cell_type1=SOURCE,
    cell_type2=TARGETS,
    means=cpdb_results["means"],
    pvals=cpdb_results["pvalues"],
    celltype_key=CELLTYPE_KEY,
    genes=GENES,
    interaction_scores=cpdb_results["interaction_scores"],
    figsize=(10, 5),
)
fair.save_fig(dirs, f"dotplot_{SOURCE}_to_{TARGETS.replace('|', '_')}".lower(), fig=dotplot)

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
