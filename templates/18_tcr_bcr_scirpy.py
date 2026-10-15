# 18 Scirpy：TCR / BCR 免疫受体库分析
#
# 【这份代码做什么】
#   把 10x VDJ 测得的 T/B 细胞受体序列与转录组数据合并，定义克隆型，
#   分析哪些 T/B 细胞发生了克隆扩增、各细胞类型的克隆多样性、V(D)J 基因使用、样本间共享的克隆，以及相似克隆组成的网络。
#   流程：读取 VDJ → 与转录组合成 MuData → 链配对质控 → 去 multichain → 定义克隆型和扩增 → 多样性、V(D)J、CDR3、共享克隆 → TCRdist 网络
# 【在流程中的位置】上游 05（+ Cell Ranger vdj 输出）→ 本步（终点分析）
# 【运行环境】sc-multi（课程 ps23）
# 【输入数据】
#   - results/05_annotate/objects/05_annotated.h5ad（含 T 或 B 细胞）
#   - 每个样本 Cell Ranger vdj 的 filtered_contig_annotations.csv（路径写在 VDJ_FILES，样本名与 01 一致）
# 【输出结果】都在 results/18_scirpy/
#   - objects/18_scirpy.h5mu
#   - tables/：18_chain_pairing_sample.csv、18_clonotypes_per_cell.csv、18_clone_sizes.csv、18_repertoire_overlap.csv
#   - figures/：18_chain_pairing_sample、18_clonal_expansion_cell_type、18_alpha_diversity_cell_type、18_umap_clonal_expansion、18_vdj_usage、18_spectratype_cell_type、18_repertoire_overlap_sample、18_clonotype_network_tcrdist
# 【怎么解读】
#   - chain_pairing：multichain 多为双细胞，已过滤；orphan = 只测到一条链；single pair 是理想情况
#   - clonal_expansion：克隆大小分级，扩增克隆集中的细胞群通常是效应或组织驻留细胞
#   - alpha_diversity（归一化 Shannon 熵）：低 = 少数克隆占主导（克隆扩增）
#   - repertoire_overlap：样本或细胞群之间共享的克隆，提示共同的抗原反应
#   - TCRdist 网络：氨基酸序列相似的克隆聚成簇，可能识别相同抗原

# %%
# 步骤：导入包，设置参数和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/18_scirpy/）
import pandas as pd
import anndata as ad
import scanpy as sc
import muon as mu
import scirpy as ir
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
VDJ_FILES = {
    "S1": "data/S1/filtered_contig_annotations.csv",
    "S2": "data/S2/filtered_contig_annotations.csv",
}
TCRDIST_CUTOFF = 15
dirs = fair.setup_step("18_scirpy")

# %%
# 步骤：读取转录组和各样本的 VDJ，细胞名格式和 05 保持一致（<barcode>-<样本>）
# 输入：IN_FILE、VDJ_FILES
# 输出：adata、airr（合并后的受体数据）
adata = sc.read_h5ad(IN_FILE)
airrs = {}
for name, path in VDJ_FILES.items():
    airrs[name] = ir.io.read_10x_vdj(path)
airr = ad.concat(airrs, index_unique="-")
airr

# %%
# 步骤：组成 MuData
# 输入：adata、airr
# 输出：mdata（模态 gex、airr）
mdata = mu.MuData({"gex": adata, "airr": airr})
mdata

# %%
# 步骤：链索引和链配对质控（图 + 表）
# 输入：mdata
# 输出：figures/18_chain_pairing_sample.pdf/.png；tables/18_chain_pairing_sample.csv
ir.pp.index_chains(mdata)
ir.tl.chain_qc(mdata)
ir.pl.group_abundance(mdata, groupby="airr:chain_pairing", target_col="gex:sample")
fair.save_fig(dirs, "chain_pairing_sample")
pd.crosstab(mdata.obs["airr:chain_pairing"], mdata.obs["gex:sample"]).to_csv(fair.out(dirs, "tables", "chain_pairing_sample.csv"))

# %%
# 步骤：去掉 multichain 细胞
# 输入：mdata
# 输出：mdata（过滤后）
mu.pp.filter_obs(mdata, "airr:chain_pairing", lambda x: x != "multichain")

# %%
# 步骤：按 CDR3 核苷酸序列完全一致定义克隆型，计算克隆扩增
# 输入：mdata
# 输出：obs['airr:clone_id']、['airr:clonal_expansion']
ir.pp.ir_dist(mdata)
ir.tl.define_clonotypes(mdata, receptor_arms="all", dual_ir="primary_only")
ir.tl.clonal_expansion(mdata)

# %%
# 步骤：保存每个细胞的克隆信息和每个克隆的大小
# 输入：mdata.obs
# 输出：tables/18_clonotypes_per_cell.csv、18_clone_sizes.csv
clone_cols = ["gex:sample", "gex:cell_type", "airr:clone_id", "airr:clone_id_size", "airr:clonal_expansion"]
mdata.obs[clone_cols].to_csv(fair.out(dirs, "tables", "clonotypes_per_cell.csv"))
mdata.obs.groupby(["airr:clone_id", "gex:sample"]).size().rename("n_cells").to_csv(fair.out(dirs, "tables", "clone_sizes.csv"))

# %%
# 步骤：各细胞类型的克隆扩增比例
# 输入：mdata
# 输出：figures/18_clonal_expansion_cell_type.pdf/.png
ir.pl.clonal_expansion(mdata, target_col="clone_id", groupby="gex:cell_type", breakpoints=(1, 2, 5))
fair.save_fig(dirs, "clonal_expansion_cell_type")

# %%
# 步骤：各细胞类型的克隆多样性
# 输入：mdata
# 输出：figures/18_alpha_diversity_cell_type.pdf/.png
ir.pl.alpha_diversity(mdata, metric="normalized_shannon_entropy", groupby="gex:cell_type")
fair.save_fig(dirs, "alpha_diversity_cell_type")

# %%
# 步骤：克隆扩增画在 UMAP 上
# 输入：mdata
# 输出：figures/18_umap_clonal_expansion.pdf/.png
mu.pl.embedding(mdata, basis="gex:umap", color="airr:clonal_expansion", show=False)
fair.save_fig(dirs, "umap_clonal_expansion")

# %%
# 步骤：V(D)J 基因组合使用
# 输入：mdata
# 输出：figures/18_vdj_usage.pdf/.png
ir.pl.vdj_usage(mdata, full_combination=False, max_ribbons=30)
fair.save_fig(dirs, "vdj_usage")

# %%
# 步骤：CDR3 长度分布
# 输入：mdata
# 输出：figures/18_spectratype_cell_type.pdf/.png
ir.pl.spectratype(mdata, color="gex:cell_type", viztype="bar")
fair.save_fig(dirs, "spectratype_cell_type")

# %%
# 步骤：样本之间克隆共享（表 + 热图）
# 输入：mdata
# 输出：tables/18_repertoire_overlap.csv；figures/18_repertoire_overlap_sample.pdf/.png
overlap, _, _ = ir.tl.repertoire_overlap(mdata, "gex:sample", inplace=False)
overlap.to_csv(fair.out(dirs, "tables", "repertoire_overlap.csv"))
grid = ir.pl.repertoire_overlap(mdata, "gex:sample", yticklabels=True, xticklabels=True)
fair.save_fig(dirs, "repertoire_overlap_sample", fig=grid)

# %%
# 步骤：氨基酸相似（TCRdist）的克隆型簇和网络图
# 输入：mdata
# 输出：obs['airr:cc_aa_tcrdist']；figures/18_clonotype_network_tcrdist.pdf/.png
ir.pp.ir_dist(mdata, metric="tcrdist", sequence="aa", cutoff=TCRDIST_CUTOFF)
ir.tl.define_clonotype_clusters(mdata, sequence="aa", metric="tcrdist", receptor_arms="all", dual_ir="any")
ir.tl.clonotype_network(mdata, min_cells=3, sequence="aa", metric="tcrdist")
ir.pl.clonotype_network(mdata, color="gex:cell_type", panel_size=(7, 7))
fair.save_fig(dirs, "clonotype_network_tcrdist")

# %%
# 步骤：保存
# 输入：mdata
# 输出：objects/18_scirpy.h5mu
mdata.write(fair.out(dirs, "objects", "scirpy.h5mu"))

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
