# 23 Metacells：把相似细胞聚合成元细胞
#
# 【这份代码做什么】
#   用 Tanay 实验室的 metacells 算法，把转录状态几乎相同的几十个细胞合并成一个“元细胞”。
#   合并后噪声更低、数据更小，适合做基因相关性、GRN 等对稀疏性敏感的分析，也可以导入 MCView 交互查看。
#   流程：UMI 分布 → 排除噪声基因和低质量细胞 → 标记 lateral 基因（细胞周期等）→ 计算元细胞 → 传递细胞类型注释 → UMAP → 保存
# 【在流程中的位置】上游 05 → 本步（终点分析，或作为其他分析的输入）
# 【运行环境】sc-metacell（课程 ps29）
# 【输入数据】
#   - results/05_annotate/objects/05_annotated.h5ad（layers['counts']；用全部基因）
# 【输出结果】都在 results/23_metacells/
#   - objects/：23_cells.h5ad（每个细胞属于哪个元细胞 obs['metacell_name']）、23_metacells.h5ad（每行一个元细胞）
#   - tables/：23_cell_to_metacell.csv、23_metacell_summary.csv、23_lateral_gene_modules.csv
#   - figures/：23_total_umis_per_cell、23_metacells_umap_cell_type
# 【怎么解读】
#   - total_umis_per_cell：两条阈值线之间的细胞被保留
#   - 每个元细胞一般包含几十个细胞；outlier 细胞不进入任何元细胞
#   - lateral 基因（细胞周期、核糖体、应激）不参与分组但仍用于离群检测；看 lateral_gene_modules 表，必要时把整个模块加入 lateral
#   - 元细胞 UMAP 按细胞类型上色：同一细胞类型的元细胞应聚在一起
# 【注意】
#   - 课程演示抽了 5 万个细胞，正式分析用全部细胞；内存不够时调小并行数

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/23_metacells/）
import matplotlib.pyplot as plt
import scanpy as sc
import metacells as mc
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
MIN_UMIS = 800                 # 课程：800–20000 UMI
MAX_UMIS = 20000
MAX_EXCLUDED_FRACTION = 0.25   # 被排除基因（线粒体等）UMI 占比上限
EXCLUDED_GENE_NAMES = ["XIST", "MALAT1"]
EXCLUDED_GENE_PATTERNS = ["MT-.*"]              # 小鼠 "mt-.*"
LATERAL_GENE_NAMES = ["AURKA", "MCM3", "MCM4", "MCM7", "MKI67", "PCNA", "RRM2", "SMC4", "TPX2",
                      "FOS", "HSP90AB1", "TXN", "DUSP1"]
LATERAL_GENE_PATTERNS = ["RP[LS].*"]
SEED = 123456
dirs = fair.setup_step("23_metacells")

# %%
# 步骤：读取数据，只保留原始 counts（metacells 要求 float32 UMI）
# 输入：IN_FILE
# 输出：full
annotated = sc.read_h5ad(IN_FILE)
full = sc.AnnData(annotated.layers["counts"].astype("float32"), obs=annotated.obs, var=annotated.var[[]])
mc.ut.set_name(full, "full")
mc.ut.top_level(full)
full

# %%
# 步骤：每个细胞 UMI 总数分布和阈值线
# 输入：full
# 输出：figures/23_total_umis_per_cell.pdf/.png
plt.hist(full.X.sum(axis=1).A1, bins=100, log=True)
plt.axvline(MIN_UMIS, color="green")
plt.axvline(MAX_UMIS, color="red")
plt.xlabel("UMIs per cell")
fair.save_fig(dirs, "total_umis_per_cell")

# %%
# 步骤：排除指定基因，并让 metacells 自动排除噪声基因
# 输入：full
# 输出：full.var['excluded_gene']
mc.pl.exclude_genes(full, excluded_gene_names=EXCLUDED_GENE_NAMES,
                    excluded_gene_patterns=EXCLUDED_GENE_PATTERNS, random_seed=SEED)

# %%
# 步骤：按 UMI 总量和被排除基因占比标记要排除的细胞
# 输入：full
# 输出：full.obs['excluded_cell']
mc.tl.compute_excluded_gene_umis(full)
mc.pl.exclude_cells(full, properly_sampled_min_cell_total=MIN_UMIS,
                    properly_sampled_max_cell_total=MAX_UMIS,
                    properly_sampled_max_excluded_genes_fraction=MAX_EXCLUDED_FRACTION)

# %%
# 步骤：取出干净数据
# 输入：full
# 输出：cells
cells = mc.pl.extract_clean_data(full, name="cells")
mc.ut.top_level(cells)
cells

# %%
# 步骤：标记 lateral 基因
# 输入：cells
# 输出：cells.var['lateral_gene']
mc.pl.mark_lateral_genes(cells, lateral_gene_names=LATERAL_GENE_NAMES,
                         lateral_gene_patterns=LATERAL_GENE_PATTERNS)

# %%
# 步骤：找与 lateral 基因相关的基因模块，保存模块表供人工检查（可把整个模块加入 lateral，见课程 ps29）
# 输入：cells
# 输出：cells.var['lateral_genes_module']；tables/23_lateral_gene_modules.csv
mc.pl.relate_to_lateral_genes(cells, random_seed=SEED)
cells.var[["lateral_gene", "lateral_genes_module"]].to_csv(fair.out(dirs, "tables", "lateral_gene_modules.csv"))

# %%
# 步骤：分治流程计算元细胞（内存不够时手动调小并行数）
# 输入：cells
# 输出：cells.obs['metacell_name']
mc.pl.set_max_parallel_piles(mc.pl.guess_max_parallel_piles(cells))
with mc.ut.progress_bar():
    mc.pl.divide_and_conquer_pipeline(cells, random_seed=SEED)

# %%
# 步骤：汇总成元细胞，并把细胞类型注释传给元细胞（取多数）
# 输入：cells
# 输出：metacells（元细胞 x 基因）
metacells = mc.pl.collect_metacells(cells, name="metacells", random_seed=SEED)
mc.tl.convey_obs_to_group(adata=cells, gdata=metacells, property_name="cell_type", method=mc.ut.most_frequent)
metacells

# %%
# 步骤：保存细胞 -> 元细胞对应表和元细胞汇总
# 输入：cells、metacells
# 输出：tables/23_cell_to_metacell.csv、23_metacell_summary.csv
cells.obs[["cell_type", "metacell_name"]].to_csv(fair.out(dirs, "tables", "cell_to_metacell.csv"))
metacells.obs.to_csv(fair.out(dirs, "tables", "metacell_summary.csv"))

# %%
# 步骤：计算元细胞 UMAP 等可视化数据（MCView 也读这些），按细胞类型上色
# 输入：cells、metacells
# 输出：metacells.obs['x']、['y']；figures/23_metacells_umap_cell_type.pdf/.png
mc.pl.compute_for_mcview(adata=cells, gdata=metacells, random_seed=SEED)
sc.pl.scatter(metacells, x="x", y="y", color="cell_type", show=False)
fair.save_fig(dirs, "metacells_umap_cell_type")

# %%
# 步骤：保存
# 输入：cells、metacells
# 输出：objects/23_cells.h5ad、23_metacells.h5ad
cells.write_h5ad(fair.out(dirs, "objects", "cells.h5ad"))
metacells.write_h5ad(fair.out(dirs, "objects", "metacells.h5ad"))

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
