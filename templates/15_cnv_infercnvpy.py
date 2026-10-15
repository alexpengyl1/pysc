# 15 inferCNVpy：推断拷贝数变异，区分恶性 / 正常细胞
#
# 【这份代码做什么】
#   按基因在染色体上的位置，把相邻基因的表达做平滑，与参考（正常）细胞比较，
#   推断每个细胞哪些染色体片段扩增或缺失。大片段 CNV 常见于肿瘤细胞，可以用来区分恶性和正常细胞。
#   流程：从 GTF 加基因坐标 → 以正常细胞为参考推断 CNV → 染色体热图 → CNV 空间聚类和打分 → 标记 tumor 簇 → 保存
# 【在流程中的位置】上游 05 → 本步（终点分析）
# 【运行环境】sc-core（课程 ps20）
# 【输入数据】
#   - results/05_annotate/objects/05_annotated.h5ad（X = log1p；obs['cell_type']）
#   - 与比对参考基因组一致的 GTF（给基因加染色体坐标）
#   - 顶部 REFERENCE_TYPES：确定是正常的细胞类型（如免疫细胞、成纤维细胞）
# 【输出结果】都在 results/15_infercnv/
#   - objects/15_infercnv.h5ad：obsm['X_cnv']、obs['cnv_score']、['cnv_leiden']、['cnv_status']
#   - tables/：15_cnv_score_per_cell.csv、15_cnv_status_by_cell_type.csv
#   - figures/：15_chromosome_heatmap_cell_type、15_chromosome_heatmap_cnv_leiden、15_cnv_umap、15_umap_cnv_status
# 【怎么解读】
#   - 染色体热图：参考细胞应基本平坦；红 = 扩增，蓝 = 缺失
#   - 成段的大范围扩增/缺失集中在某些上皮或肿瘤群 → 判为恶性；cnv_score 高的 cnv_leiden 群
#   - cnv_status_by_cell_type：每种细胞类型里被判为 tumor 的细胞数
# 【注意】
#   - TUMOR_CLUSTERS 要看完热图后人工填写
#   - 非肿瘤样本（如放射损伤组织）一般不会有大片段 CNV，结果只能作为排除恶变的辅助证据

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/15_infercnv/）
import pandas as pd
import scanpy as sc
import infercnvpy as cnv
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
GTF = "data/genes.gtf"
REFERENCE_TYPES = ["T", "B", "Myeloid", "Endothelial"]   # 当作正常参考的细胞类型
WINDOW_SIZE = 250
TUMOR_CLUSTERS = ["1", "4"]                               # 看完 cnv_leiden 热图后填写
dirs = fair.setup_step("15_infercnv")

# %%
# 步骤：读取数据，从 GTF 加基因组坐标
# 输入：IN_FILE、GTF
# 输出：adata.var 新增 chromosome、start、end
adata = sc.read_h5ad(IN_FILE)
cnv.io.genomic_position_from_gtf(GTF, adata, gtf_gene_id="gene_name")
adata.var[["chromosome", "start", "end"]].head()

# %%
# 步骤：以正常细胞为参考推断 CNV
# 输入：adata
# 输出：adata.obsm['X_cnv']
cnv.tl.infercnv(adata, reference_key="cell_type", reference_cat=REFERENCE_TYPES, window_size=WINDOW_SIZE)

# %%
# 步骤：按细胞类型画染色体热图
# 输入：adata.obsm['X_cnv']
# 输出：figures/15_chromosome_heatmap_cell_type.pdf/.png
cnv.pl.chromosome_heatmap(adata, groupby="cell_type", show=False)
fair.save_fig(dirs, "chromosome_heatmap_cell_type")

# %%
# 步骤：在 CNV 空间聚类，计算每个细胞的 CNV 分数
# 输入：adata.obsm['X_cnv']
# 输出：obs['cnv_leiden']、obs['cnv_score']、obsm['X_cnv_umap']
cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata)
cnv.tl.umap(adata)
cnv.tl.cnv_score(adata)

# %%
# 步骤：按 CNV 簇画染色体热图
# 输入：adata
# 输出：figures/15_chromosome_heatmap_cnv_leiden.pdf/.png
cnv.pl.chromosome_heatmap(adata, groupby="cnv_leiden", dendrogram=True, show=False)
fair.save_fig(dirs, "chromosome_heatmap_cnv_leiden")

# %%
# 步骤：CNV 空间的 UMAP（CNV 簇、分数、细胞类型）
# 输入：adata
# 输出：figures/15_cnv_umap.pdf/.png
cnv.pl.umap(adata, color=["cnv_leiden", "cnv_score", "cell_type"], show=False)
fair.save_fig(dirs, "cnv_umap")

# %%
# 步骤：把 CNV 明显的簇标为 tumor
# 输入：adata.obs['cnv_leiden']、TUMOR_CLUSTERS
# 输出：obs['cnv_status']；figures/15_umap_cnv_status.pdf/.png
adata.obs["cnv_status"] = "normal"
adata.obs.loc[adata.obs["cnv_leiden"].isin(TUMOR_CLUSTERS), "cnv_status"] = "tumor"
sc.pl.umap(adata, color=["cnv_status", "cell_type"], show=False)
fair.save_fig(dirs, "umap_cnv_status")

# %%
# 步骤：保存每个细胞的 CNV 结果和各细胞类型的 tumor/normal 计数
# 输入：adata.obs
# 输出：tables/15_cnv_score_per_cell.csv、15_cnv_status_by_cell_type.csv
adata.obs[["cell_type", "cnv_leiden", "cnv_score", "cnv_status"]].to_csv(fair.out(dirs, "tables", "cnv_score_per_cell.csv"))
pd.crosstab(adata.obs["cell_type"], adata.obs["cnv_status"]).to_csv(fair.out(dirs, "tables", "cnv_status_by_cell_type.csv"))

# %%
# 步骤：保存
# 输入：adata
# 输出：objects/15_infercnv.h5ad
adata.write(fair.out(dirs, "objects", "infercnv.h5ad"))

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
