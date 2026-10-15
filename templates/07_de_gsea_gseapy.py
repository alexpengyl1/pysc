# 07 某个细胞类型内：条件间差异基因 + GSEA + Enrichr 富集分析
#
# 【这份代码做什么】
#   在一个指定的细胞类型里比较两个条件（CASE vs CONTROL），找差异基因，
#   再用两种方法看哪些通路改变：GSEA prerank（用全部基因的排序）和 Enrichr（只用显著的上调/下调基因）。
#   流程：下载基因集库 → 取出细胞类型 → Wilcoxon 差异分析 → GSEA prerank → Enrichr（上调、下调分开）→ 作图
# 【在流程中的位置】上游 05 → 本步（终点分析；换 CELLTYPE 可对每个细胞类型重复运行）
# 【运行环境】gseapy（课程 ps12）；需要联网下载基因集
# 【输入数据】
#   - results/05_annotate/objects/05_annotated.h5ad（X = log1p；obs['cell_type']、obs['condition']）
#   - 顶部 CELLTYPE、CASE、CONTROL、基因集库名（GO / KEGG / MSigDB）
# 【输出结果】都在 results/07_gseapy/，文件名带细胞类型（如 07_deg_cd14_mono.csv）
#   - tables/：07_deg_<细胞类型>.csv（全部基因的差异结果）、07_prerank_<细胞类型>.csv、07_enrichr_up_<细胞类型>.csv、07_enrichr_down_<细胞类型>.csv
#   - figures/：07_deg_top25_<细胞类型>、07_prerank_top5_<细胞类型>、07_enrichr_up_dotplot_<细胞类型>、07_enrichr_down_dotplot_<细胞类型>
# 【怎么解读】
#   - 差异基因表：scores > 0 = 在 CASE 组上调；看 pvals_adj 和 logfoldchanges
#   - prerank：NES > 0 = 通路在 CASE 组整体上调，NES < 0 = 下调；FDR q-val < 0.25 是 GSEA 常用阈值
#   - 富集曲线：峰在左侧 = 通路基因集中在上调端
#   - Enrichr 点图：Adjusted P-value < 0.05 的通路；上调和下调分开解读
# 【注意】
#   - 逐细胞的检验把细胞当作重复，p 值会偏小；有多个样本时更稳妥的是按样本做 pseudobulk

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/07_gseapy/）
import scanpy as sc
import gseapy as gp
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
CELLTYPE = "CD14 Mono"
CONDITION_KEY = "condition"
CASE, CONTROL = "STIM", "CTRL"
ORGANISM = "Human"                         # 或 "Mouse"
LIBRARY = "GO_Biological_Process_2023"      # gp.get_library_name(organism=ORGANISM) 查看全部
SEED = 0
TAG = CELLTYPE.lower().replace(" ", "_")   # 文件名里的细胞类型标签
dirs = fair.setup_step("07_gseapy")

# %%
# 步骤：从 Enrichr 下载基因集库
# 输入：LIBRARY、ORGANISM
# 输出：gene_sets（字典：通路名 -> 基因列表）
gene_sets = gp.get_library(name=LIBRARY, organism=ORGANISM)
len(gene_sets)

# %%
# 步骤：读取数据，取出一个细胞类型
# 输入：IN_FILE
# 输出：bdata（该细胞类型的子集）；打印两组细胞数
adata = sc.read_h5ad(IN_FILE)
bdata = adata[adata.obs["cell_type"] == CELLTYPE].copy()
bdata.obs[CONDITION_KEY].value_counts()

# %%
# 步骤：CASE vs CONTROL 差异分析（Wilcoxon）
# 输入：bdata
# 输出：deg（差异表）；tables/07_deg_<细胞类型>.csv
sc.tl.rank_genes_groups(bdata, groupby=CONDITION_KEY, groups=[CASE], reference=CONTROL,
                        method="wilcoxon", use_raw=False)
deg = sc.get.rank_genes_groups_df(bdata, group=CASE)
deg.to_csv(fair.out(dirs, "tables", f"deg_{TAG}.csv"), index=False)
deg.head()

# %%
# 步骤：前 25 个差异基因图
# 输入：bdata.uns['rank_genes_groups']
# 输出：figures/07_deg_top25_<细胞类型>.pdf/.png
sc.pl.rank_genes_groups(bdata, n_genes=25, sharey=False, show=False)
fair.save_fig(dirs, f"deg_top25_{TAG}")

# %%
# 步骤：按 Wilcoxon 分数排序做 GSEA prerank
# 输入：deg（names, scores）
# 输出：prerank（结果对象）；tables/07_prerank_<细胞类型>.csv
prerank = gp.prerank(rnk=deg[["names", "scores"]], gene_sets=gene_sets,
                     permutation_num=1000, seed=SEED, outdir=None)
prerank.res2d.to_csv(fair.out(dirs, "tables", f"prerank_{TAG}.csv"), index=False)
prerank.res2d.head(10)

# %%
# 步骤：前 5 个通路的富集曲线
# 输入：prerank
# 输出：figures/07_prerank_top5_<细胞类型>.pdf/.png
prerank.plot(terms=prerank.res2d.Term[:5])
fair.save_fig(dirs, f"prerank_top5_{TAG}")

# %%
# 步骤：显著上调 / 下调基因分别做 Enrichr（背景 = 检测到的全部基因）
# 输入：deg
# 输出：enr_up、enr_down；tables/07_enrichr_up_<细胞类型>.csv、07_enrichr_down_<细胞类型>.csv
sig = deg[deg.pvals_adj < 0.05]
background = bdata.var_names.tolist()
enr_up = gp.enrichr(gene_list=sig[sig.logfoldchanges > 0].names.tolist(), gene_sets=gene_sets,
                    background=background, outdir=None)
enr_down = gp.enrichr(gene_list=sig[sig.logfoldchanges < 0].names.tolist(), gene_sets=gene_sets,
                      background=background, outdir=None)
enr_up.res2d.to_csv(fair.out(dirs, "tables", f"enrichr_up_{TAG}.csv"), index=False)
enr_down.res2d.to_csv(fair.out(dirs, "tables", f"enrichr_down_{TAG}.csv"), index=False)

# %%
# 步骤：上调基因 Enrichr 点图
# 输入：enr_up
# 输出：figures/07_enrichr_up_dotplot_<细胞类型>.pdf/.png
gp.dotplot(enr_up.res2d, title=f"{CELLTYPE} up in {CASE}", figsize=(4, 6))
fair.save_fig(dirs, f"enrichr_up_dotplot_{TAG}")

# %%
# 步骤：下调基因 Enrichr 点图
# 输入：enr_down
# 输出：figures/07_enrichr_down_dotplot_<细胞类型>.pdf/.png
gp.dotplot(enr_down.res2d, title=f"{CELLTYPE} down in {CASE}", figsize=(4, 6))
fair.save_fig(dirs, f"enrichr_down_dotplot_{TAG}")

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
