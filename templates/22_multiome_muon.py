# 22 muon：10x Multiome（RNA + ATAC）分析与整合
#
# 【这份代码做什么】
#   同一个细胞同时测了基因表达（RNA）和染色质开放性（ATAC）时，先分别对两个模态做质控、标准化、降维和聚类，
#   再用 MOFA+（找共享和特有的变异因子）和 WNN（按细胞加权组合两个模态的邻居图）整合成一个联合表示。
#   流程：RNA 质控/聚类/marker → ATAC 质控（fragment 长度、核小体信号、TSS 富集）/聚类 → 同步细胞 → MOFA+ → WNN UMAP
# 【在流程中的位置】独立分析（直接读 Cell Ranger ARC 输出，不依赖 01–05）
# 【运行环境】sc-multi（课程 ps28）
# 【输入数据】
#   - Cell Ranger ARC 输出的 filtered_feature_bc_matrix.h5，同一文件夹里放 atac_fragments.tsv.gz 和 .tbi 索引
# 【输出结果】都在 results/22_muon/
#   - objects/22_multiome.h5mu：rna、atac 两个模态；obsm['X_mofa']、WNN UMAP；objects/22_mofa_model.hdf5
#   - tables/：22_qc_cells_kept.csv、22_markers_rna_leiden.csv
#   - figures/：22_rna_qc_violin、22_rna_umap_leiden、22_atac_qc_histogram、22_atac_fragment_length、22_atac_nucleosome_signal、22_atac_tss_enrichment、22_atac_umap_leiden、22_wnn_umap
# 【怎么解读】
#   - RNA 质控：同 02
#   - fragment 长度图：应能看到 < 147 bp 的无核小体峰和约 200 bp 的单核小体峰（周期性）
#   - nucleosome_signal：过高的细胞质量差，去掉
#   - TSS 富集：转录起始位点处的峰越尖锐，数据质量越好
#   - MOFA 因子：既能捕获两个模态共享的变异，也能捕获模态特有的变异
#   - WNN UMAP：分别按 RNA 聚类和 ATAC 聚类上色，两者一致说明整合合理
# 【注意】
#   - 小鼠线粒体前缀用 mt-

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/22_muon/）
import pandas as pd
import scanpy as sc
import muon as mu
from muon import atac as ac
import fair_utils as fair

H5_FILE = "data/multiome/filtered_feature_bc_matrix.h5"
MT_PREFIX = "MT-"
RNA_MIN_GENES, RNA_MAX_GENES, RNA_MAX_PCT_MT = 200, 8000, 10
ATAC_MIN_COUNTS, ATAC_MAX_COUNTS = 1000, 80000
N_FACTORS = 20
RESOLUTION = 0.5
SEED = 0
dirs = fair.setup_step("22_muon")

# %%
# 步骤：读取多组学数据
# 输入：H5_FILE（fragments 文件自动关联）
# 输出：mdata（模态 rna、atac）
mdata = mu.read_10x_h5(H5_FILE)
mdata.var_names_make_unique()
cells_before = {"rna": mdata["rna"].n_obs, "atac": mdata["atac"].n_obs}
mdata

# %%
# 步骤：RNA 质控指标和小提琴图
# 输入：mdata['rna']
# 输出：rna；figures/22_rna_qc_violin.pdf/.png
rna = mdata.mod["rna"]
rna.var["mt"] = rna.var_names.str.startswith(MT_PREFIX)
sc.pp.calculate_qc_metrics(rna, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(rna, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], multi_panel=True, show=False)
fair.save_fig(dirs, "rna_qc_violin")

# %%
# 步骤：RNA 细胞过滤
# 输入：rna
# 输出：rna（过滤后）
mu.pp.filter_obs(rna, "n_genes_by_counts", lambda x: (x >= RNA_MIN_GENES) & (x < RNA_MAX_GENES))
mu.pp.filter_obs(rna, "pct_counts_mt", lambda x: x < RNA_MAX_PCT_MT)

# %%
# 步骤：RNA 标准化、高变基因、PCA、聚类、UMAP
# 输入：rna
# 输出：rna.obs['leiden']、obsm['X_pca']、['X_umap']；figures/22_rna_umap_leiden.pdf/.png
rna.layers["counts"] = rna.X.copy()
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
rna.layers["log1p"] = rna.X.copy()
sc.pp.highly_variable_genes(rna, min_mean=0.02, max_mean=4, min_disp=0.5)
sc.pp.scale(rna, max_value=10)
sc.tl.pca(rna, random_state=SEED)
sc.pp.neighbors(rna, n_pcs=20, random_state=SEED)
sc.tl.leiden(rna, resolution=RESOLUTION, flavor="igraph", n_iterations=2, directed=False, random_state=SEED)
sc.tl.umap(rna, random_state=SEED)
sc.pl.umap(rna, color="leiden", legend_loc="on data", show=False)
fair.save_fig(dirs, "rna_umap_leiden")

# %%
# 步骤：RNA 各群 marker（用 log1p 层）
# 输入：rna
# 输出：tables/22_markers_rna_leiden.csv
sc.tl.rank_genes_groups(rna, "leiden", method="wilcoxon", layer="log1p", use_raw=False)
sc.get.rank_genes_groups_df(rna, group=None).to_csv(fair.out(dirs, "tables", "markers_rna_leiden.csv"), index=False)

# %%
# 步骤：ATAC 质控指标分布
# 输入：mdata['atac']
# 输出：atac；figures/22_atac_qc_histogram.pdf/.png
atac = mdata.mod["atac"]
sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)
mu.pl.histogram(atac, ["n_genes_by_counts", "total_counts"])
fair.save_fig(dirs, "atac_qc_histogram")

# %%
# 步骤：ATAC 峰和细胞过滤
# 输入：atac
# 输出：atac（过滤后）
mu.pp.filter_var(atac, "n_cells_by_counts", lambda x: x >= 10)
mu.pp.filter_obs(atac, "total_counts", lambda x: (x >= ATAC_MIN_COUNTS) & (x <= ATAC_MAX_COUNTS))

# %%
# 步骤：fragment 长度分布
# 输入：atac（fragments 文件）
# 输出：figures/22_atac_fragment_length.pdf/.png
ac.pl.fragment_histogram(atac, region="chr1:1-2000000")
fair.save_fig(dirs, "atac_fragment_length")

# %%
# 步骤：核小体信号
# 输入：atac
# 输出：atac.obs['nucleosome_signal']；figures/22_atac_nucleosome_signal.pdf/.png
ac.tl.nucleosome_signal(atac, n=1e6)
mu.pl.histogram(atac, "nucleosome_signal")
fair.save_fig(dirs, "atac_nucleosome_signal")

# %%
# 步骤：TSS 富集
# 输入：mdata
# 输出：atac.obs['tss_score']；figures/22_atac_tss_enrichment.pdf/.png
tss = ac.tl.tss_enrichment(mdata, n_tss=1000)
ac.pl.tss_enrichment(tss)
fair.save_fig(dirs, "atac_tss_enrichment")

# %%
# 步骤：ATAC 标准化、高变峰、PCA、聚类、UMAP
# 输入：atac
# 输出：atac.obs['leiden']、obsm['X_pca']、['X_umap']；figures/22_atac_umap_leiden.pdf/.png
atac.layers["counts"] = atac.X.copy()
sc.pp.normalize_total(atac, target_sum=1e4)
sc.pp.log1p(atac)
sc.pp.highly_variable_genes(atac, min_mean=0.05, max_mean=1.5, min_disp=0.5)
sc.pp.scale(atac, max_value=10)
sc.tl.pca(atac, random_state=SEED)
sc.pp.neighbors(atac, n_pcs=20, random_state=SEED)
sc.tl.leiden(atac, resolution=RESOLUTION, flavor="igraph", n_iterations=2, directed=False, random_state=SEED)
sc.tl.umap(atac, random_state=SEED)
sc.pl.umap(atac, color="leiden", legend_loc="on data", show=False)
fair.save_fig(dirs, "atac_umap_leiden")

# %%
# 步骤：同步两个模态，只保留两边都有的细胞，记录细胞数
# 输入：mdata
# 输出：mdata（细胞对齐）；tables/22_qc_cells_kept.csv
mdata.update()
mu.pp.intersect_obs(mdata)
pd.DataFrame({"before": cells_before, "after": {"rna": mdata["rna"].n_obs, "atac": mdata["atac"].n_obs}}).to_csv(
    fair.out(dirs, "tables", "qc_cells_kept.csv"))

# %%
# 步骤：MOFA+ 因子整合
# 输入：mdata
# 输出：mdata.obsm['X_mofa']；objects/22_mofa_model.hdf5
mu.tl.mofa(mdata, n_factors=N_FACTORS, outfile=str(fair.out(dirs, "objects", "mofa_model.hdf5")))

# %%
# 步骤：WNN 加权邻居图和 UMAP（两个模态先各自重算邻居图）
# 输入：mdata
# 输出：mdata.obsm['X_umap']；figures/22_wnn_umap.pdf/.png
sc.pp.neighbors(mdata["rna"])
sc.pp.neighbors(mdata["atac"])
mu.pp.neighbors(mdata, key_added="wnn")
mu.tl.umap(mdata, neighbors_key="wnn", random_state=SEED)
mu.pl.umap(mdata, color=["rna:leiden", "atac:leiden"], show=False)
fair.save_fig(dirs, "wnn_umap")

# %%
# 步骤：保存
# 输入：mdata
# 输出：objects/22_multiome.h5mu
mdata.write(fair.out(dirs, "objects", "multiome.h5mu"))

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
