# 20 scDRS：找与疾病 / 性状相关的细胞
#
# 【这份代码做什么】
#   把 GWAS 得到的疾病（或性状）相关基因集映射到单细胞数据上，给每个细胞打一个疾病相关得分，
#   并与大量表达量匹配的随机基因集比较得到显著性，再判断哪些细胞类型与疾病相关、类型内部是否存在异质性。
#   流程：标准化 → 协变量 → 读基因集 → 回归协变量和基因分箱 → 逐性状打分 → UMAP → 细胞类型水平检验
# 【在流程中的位置】上游 05 → 本步（终点分析）
# 【运行环境】sc-core（课程 ps26）
# 【输入数据】
#   - results/05_annotate/objects/05_annotated.h5ad（layers['counts']）
#   - 基因集文件 .gs：两列 TRAIT, GENESET，基因格式 "GENE1:权重,GENE2:权重"；
#     由 GWAS summary statistics 跑 MAGMA 后用 scdrs munge-gs 生成，或下载 scDRS 提供的 74 个性状
# 【输出结果】都在 results/20_scdrs/，文件名带性状名
#   - tables/：20_score_<性状>.csv（每个细胞的得分、p 值和 1000 组对照得分）、20_group_stats_<性状>.csv
#   - figures/20_umap_scores
# 【怎么解读】
#   - norm_score：每个细胞的标准化疾病得分，越高越相关；UMAP 上集中在哪些细胞
#   - group_stats 表：assoc_mcp 小 = 该细胞类型与疾病显著相关；hetero_mcp 小 = 类型内部细胞间的关联有异质性（可能是某个亚群）
#   - 用 Height 等无关性状作阴性对照，确认结果不是普遍存在的
# 【注意】
#   - 小鼠数据用 load_gs 的同源转换参数

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/20_scdrs/）
import pandas as pd
import scanpy as sc
import scdrs
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
GS_FILE = "data/geneset.gs"
SPECIES = "human"          # 表达数据物种：human / mouse
GS_SPECIES = "human"       # 基因集物种（MAGMA 结果一般是 human）
N_CTRL = 1000
dirs = fair.setup_step("20_scdrs")

# %%
# 步骤：读取数据，X 换成原始 counts，做 scDRS 推荐的过滤和标准化
# 输入：IN_FILE
# 输出：adata（X = log 标准化后的 counts）
adata = sc.read_h5ad(IN_FILE)
adata.X = adata.layers["counts"].copy()
sc.pp.filter_cells(adata, min_genes=250)
sc.pp.filter_genes(adata, min_cells=50)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# %%
# 步骤：协变量表（常数项 + 每细胞检测到的基因数）
# 输入：adata.obs
# 输出：cov（细胞 x 协变量）
cov = pd.DataFrame({"const": 1, "n_genes": adata.obs["n_genes"]}, index=adata.obs_names)

# %%
# 步骤：读取基因集，必要时做同源转换，并与表达基因取交集
# 输入：GS_FILE
# 输出：gene_sets（性状 -> (基因列表, 权重列表)）
gene_sets = scdrs.util.load_gs(GS_FILE, src_species=GS_SPECIES, dst_species=SPECIES,
                               to_intersect=adata.var_names)
list(gene_sets.keys())

# %%
# 步骤：回归协变量，按表达均值/方差给基因分箱
# 输入：adata、cov
# 输出：adata（就地修改）
scdrs.preprocess(adata, cov=cov, n_mean_bin=20, n_var_bin=20, copy=False)

# %%
# 步骤：逐个性状给每个细胞打分
# 输入：adata、gene_sets
# 输出：scores（性状 -> 得分表）；每个性状一个 tables/20_score_<性状>.csv
scores = {}
for trait, (genes, weights) in gene_sets.items():
    scores[trait] = scdrs.score_cell(adata, genes, gene_weight=weights, ctrl_match_key="mean_var",
                                     n_ctrl=N_CTRL, weight_opt="vs", return_ctrl_norm_score=True)
    scores[trait].to_csv(fair.out(dirs, "tables", f"score_{trait}.csv"))

# %%
# 步骤：得分画在 UMAP 上
# 输入：scores
# 输出：figures/20_umap_scores.pdf/.png
for trait in scores:
    adata.obs[trait] = scores[trait]["norm_score"]
sc.pl.umap(adata, color=list(scores), cmap="RdBu_r", vmin=-5, vmax=5, show=False)
fair.save_fig(dirs, "umap_scores")

# %%
# 步骤：细胞类型水平的关联和异质性检验
# 输入：adata、scores
# 输出：每个性状一个 tables/20_group_stats_<性状>.csv
for trait in scores:
    stats = scdrs.method.downstream_group_analysis(adata, df_full_score=scores[trait], group_cols=["cell_type"])
    stats["cell_type"].to_csv(fair.out(dirs, "tables", f"group_stats_{trait}.csv"))

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
