# 17 drug2cell：寻找细胞类型特异的药物靶点信号
#
# 【这份代码做什么】
#   用 ChEMBL 数据库里每种药物的靶点基因集，给每个细胞打分（靶点基因的平均表达），
#   找出哪些药物的靶点在哪类细胞里特异高表达，提示药物可能作用的细胞，或潜在的脱靶效应。
#   流程：药物靶点打分 → 在打分空间找各细胞类型特异药物 → 点图 → 细胞类型 marker 与靶点集的超几何检验
# 【在流程中的位置】上游 05 → 本步（终点分析）
# 【运行环境】sc-core（课程 ps22）
# 【输入数据】
#   - results/05_annotate/objects/05_annotated.h5ad（X = log1p；obs['cell_type']）
#   - ChEMBL 靶点库：drug2cell 包自带
# 【输出结果】都在 results/17_drug2cell/
#   - objects/17_drug_scores.h5ad：细胞 x 药物 打分矩阵
#   - tables/：17_drug_markers_cell_type.csv、17_hypergeom_<细胞类型>.csv
#   - figures/17_drug_dotplot_cell_type
# 【怎么解读】
#   - 点图：某药物在某细胞类型的打分高 → 该药物可能作用于这类细胞（包括潜在脱靶）
#   - drug_markers 表：每个细胞类型特异的药物，scores 高、pvals_adj 小优先
#   - hypergeom 表：该细胞类型的 marker 基因与药物靶点集的重叠是否显著，作为另一条证据
# 【注意】
#   - 靶点库是人类基因名，小鼠数据先换成人类同源基因（见 09_ccc_cellphonedb.py 的同源转换块）
#   - 结果只是基于表达的提示，需要实验或文献支持

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/17_drug2cell/）
import scanpy as sc
import drug2cell as d2c
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
CELLTYPE_KEY = "cell_type"
N_TOP = 5
dirs = fair.setup_step("17_drug2cell")

# %%
# 步骤：读取数据
# 输入：IN_FILE
# 输出：adata
adata = sc.read_h5ad(IN_FILE)

# %%
# 步骤：用 ChEMBL 药物靶点集给每个细胞打分，保存打分矩阵
# 输入：adata.X（log1p）
# 输出：drug_scores（细胞 x 药物 的 AnnData）；objects/17_drug_scores.h5ad
d2c.score(adata, use_raw=False)
drug_scores = adata.uns.pop("drug2cell")
drug_scores.write(fair.out(dirs, "objects", "drug_scores.h5ad"))
drug_scores

# %%
# 步骤：在药物打分空间找各细胞类型特异的药物
# 输入：drug_scores
# 输出：tables/17_drug_markers_cell_type.csv
sc.tl.rank_genes_groups(drug_scores, groupby=CELLTYPE_KEY, method="wilcoxon")
sc.get.rank_genes_groups_df(drug_scores, group=None).to_csv(fair.out(dirs, "tables", "drug_markers_cell_type.csv"), index=False)

# %%
# 步骤：各细胞类型前 N_TOP 个特异药物的点图
# 输入：drug_scores
# 输出：figures/17_drug_dotplot_cell_type.pdf/.png
sc.pl.rank_genes_groups_dotplot(drug_scores, n_genes=N_TOP, swap_axes=True, dendrogram=False, show=False)
fair.save_fig(dirs, "drug_dotplot_cell_type")

# %%
# 步骤：基因表达空间的 marker，用于过表征分析
# 输入：adata
# 输出：adata.uns['rank_genes_groups']
sc.tl.rank_genes_groups(adata, groupby=CELLTYPE_KEY, method="wilcoxon", use_raw=False)

# %%
# 步骤：marker 与药物靶点集的超几何检验
# 输入：adata
# 输出：每个细胞类型一个 tables/17_hypergeom_<细胞类型>.csv
overrep = d2c.hypergeometric(adata)
for celltype, table in overrep.items():
    table.to_csv(fair.out(dirs, "tables", f"hypergeom_{celltype.replace(' ', '_')}.csv"))

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
