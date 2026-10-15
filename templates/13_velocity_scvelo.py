# 13 RNA velocity（scVelo）
#
# 【这份代码做什么】
#   利用未剪接（unspliced）和已剪接（spliced）mRNA 的比例，推断每个细胞的基因表达正在朝哪个方向变化（RNA 速度），
#   从而预测细胞状态的转变方向，并估计内部的“潜在时间”和驱动这些变化的基因。
#   流程：读入 loom 并合并 → 剪接比例 → 标准化和邻域矩 → 动力学模型 → 速度流线图 → 置信度 → latent time → 驱动基因
# 【在流程中的位置】上游 05（+ velocyto 生成的 loom）→ 本步 → 下游 14_fate_cellrank（可选）
# 【运行环境】sc-velo（课程 ps18）
# 【输入数据】
#   - 每个样本的 loom（spliced / unspliced counts）；loom 要先用 velocyto 从 Cell Ranger 输出生成（需要 BAM、基因组 GTF、可选 repeat mask GTF）：
#     pip install velocyto（需要 samtools），然后 velocyto run10x -m repeat_msk.gtf <cellranger_sample_dir> genes.gtf，
#     生成的文件在 <cellranger_sample_dir>/velocyto/<sample>.loom
#   - results/05_annotate/objects/05_annotated.h5ad
# 【输出结果】都在 results/13_scvelo/
#   - objects/13_velocity.h5ad（14 可用）
#   - tables/：13_spliced_proportions.csv、13_velocity_confidence.csv、13_top_likelihood_genes.csv、13_dynamical_genes_cell_type.csv
#   - figures/：13_proportions_cell_type、13_velocity_stream_cell_type、13_velocity_confidence_length、13_latent_time、13_heatmap_top300_latent_time、13_phase_portraits_top10
# 【怎么解读】
#   - proportions：unspliced 比例通常约 10–25%，太低时速度不可靠
#   - 流线图：箭头 = 推断的细胞状态变化方向；要结合单基因相图核对，不能只看箭头
#   - velocity_confidence：低的区域方向不可信
#   - latent_time：基于转录动力学的“内部时间”
#   - 驱动基因：fit_likelihood 高的基因动态最明显；相图（spliced vs unspliced）呈现清晰的弧形 = 拟合好
# 【注意】
#   - 终末分化、静止的细胞群速度很弱，箭头方向容易误导

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/13_scvelo/）
import anndata as ad
import scanpy as sc
import scvelo as scv
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
LOOMS = {"S1": "data/velocyto/S1.loom", "S2": "data/velocyto/S2.loom"}
MIN_SHARED_COUNTS = 20
N_TOP_GENES = 2000
MODE = "dynamical"      # 赶时间可改 "stochastic"，并删掉 recover_dynamics 和 latent_time 相关块
N_JOBS = 4
scv.settings.verbosity = 3
dirs = fair.setup_step("13_scvelo")

# %%
# 步骤：读取注释数据和 loom；loom 细胞名改成和 05 一致的格式（<barcode>-1-<样本>）
# 输入：IN_FILE、LOOMS
# 输出：adata、ldata（所有样本 loom 合并）
adata = sc.read_h5ad(IN_FILE)
looms = {}
for name, path in LOOMS.items():
    looms[name] = scv.read(path, cache=True)
    looms[name].var_names_make_unique()
    looms[name].obs_names = [bc.split(":")[1].replace("x", "-1") + "-" + name for bc in looms[name].obs_names]
ldata = ad.concat(looms)
ldata

# %%
# 步骤：把 spliced/unspliced 合并进注释数据
# 输入：adata、ldata
# 输出：adata（新增 layers['spliced']、['unspliced']）
adata = scv.utils.merge(adata, ldata)
adata

# %%
# 步骤：各细胞类型的 spliced / unspliced 比例（图 + 表）
# 输入：adata
# 输出：figures/13_proportions_cell_type.pdf/.png；tables/13_spliced_proportions.csv
scv.pl.proportions(adata, groupby="cell_type", show=False)
fair.save_fig(dirs, "proportions_cell_type")
proportions = adata.obs[["cell_type"]].assign(
    spliced=adata.layers["spliced"].sum(axis=1).A1, unspliced=adata.layers["unspliced"].sum(axis=1).A1)
proportions.groupby("cell_type")[["spliced", "unspliced"]].sum().to_csv(fair.out(dirs, "tables", "spliced_proportions.csv"))

# %%
# 步骤：过滤标准化，计算邻域矩
# 输入：adata
# 输出：adata.layers['Ms']、['Mu']
scv.pp.filter_and_normalize(adata, min_shared_counts=MIN_SHARED_COUNTS, n_top_genes=N_TOP_GENES)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# %%
# 步骤：动力学模型估计速度（耗时）
# 输入：adata
# 输出：adata.layers['velocity']、uns['velocity_graph']
scv.tl.recover_dynamics(adata, n_jobs=N_JOBS)
scv.tl.velocity(adata, mode=MODE)
scv.tl.velocity_graph(adata)

# %%
# 步骤：速度流线图
# 输入：adata（obsm['X_umap']）
# 输出：figures/13_velocity_stream_cell_type.pdf/.png
scv.pl.velocity_embedding_stream(adata, basis="umap", color="cell_type", show=False)
fair.save_fig(dirs, "velocity_stream_cell_type")

# %%
# 步骤：速度置信度和速度长度（图 + 表）
# 输入：adata
# 输出：figures/13_velocity_confidence_length.pdf/.png；tables/13_velocity_confidence.csv
scv.tl.velocity_confidence(adata)
scv.pl.scatter(adata, c=["velocity_length", "velocity_confidence"], cmap="coolwarm", perc=[5, 95], show=False)
fair.save_fig(dirs, "velocity_confidence_length")
adata.obs[["cell_type", "velocity_length", "velocity_confidence"]].to_csv(fair.out(dirs, "tables", "velocity_confidence.csv"))

# %%
# 步骤：latent time
# 输入：adata
# 输出：obs['latent_time']；figures/13_latent_time.pdf/.png
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color="latent_time", color_map="gnuplot", show=False)
fair.save_fig(dirs, "latent_time")

# %%
# 步骤：按拟合似然排序的驱动基因（表）
# 输入：adata.var['fit_likelihood']
# 输出：top_genes；tables/13_top_likelihood_genes.csv
top_genes = adata.var["fit_likelihood"].sort_values(ascending=False)
top_genes.to_csv(fair.out(dirs, "tables", "top_likelihood_genes.csv"))
top_genes = top_genes.index

# %%
# 步骤：前 300 个驱动基因沿 latent time 的热图
# 输入：adata、top_genes
# 输出：figures/13_heatmap_top300_latent_time.pdf/.png
scv.pl.heatmap(adata, var_names=top_genes[:300], sortby="latent_time", col_color="cell_type", n_convolve=100, show=False)
fair.save_fig(dirs, "heatmap_top300_latent_time")

# %%
# 步骤：前 10 个驱动基因的相图
# 输入：adata、top_genes
# 输出：figures/13_phase_portraits_top10.pdf/.png
scv.pl.scatter(adata, basis=top_genes[:10], ncols=5, color="cell_type", frameon=False, show=False)
fair.save_fig(dirs, "phase_portraits_top10")

# %%
# 步骤：各细胞类型的特异动态基因
# 输入：adata
# 输出：tables/13_dynamical_genes_cell_type.csv
scv.tl.rank_dynamical_genes(adata, groupby="cell_type")
scv.get_df(adata, "rank_dynamical_genes/names").to_csv(fair.out(dirs, "tables", "dynamical_genes_cell_type.csv"))

# %%
# 步骤：保存
# 输入：adata
# 输出：objects/13_velocity.h5ad
adata.write(fair.out(dirs, "objects", "velocity.h5ad"))

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
