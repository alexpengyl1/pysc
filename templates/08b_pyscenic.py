# 08b pySCENIC：推断转录因子调控网络（regulon）和每个细胞的活性
#
# 【这份代码做什么】
#   找出各转录因子（TF）调控哪些靶基因，并计算每个细胞里这些调控模块（regulon）的活性。分三个阶段：
#     阶段 I  GRNBoost2：根据共表达推断 TF → 靶基因 关系
#     阶段 II cisTarget：只保留靶基因启动子附近有该 TF 结合 motif 的关系，得到 regulon
#     阶段 III AUCell：计算每个细胞中每个 regulon 的活性（AUC）
# 【在流程中的位置】上游 08a → 本步 → 下游 08c_scenic_plot
# 【运行环境】pyscenic（Python 3.7；课程 ps13）
# 【输入数据】
#   - results/08_scenic/objects/08_expr_counts.csv（08a 导出）
#   - data/cistarget/ 里的数据库（https://resources.aertslab.org/cistarget/）：
#     TF 列表 allTFs_hg38.txt / allTFs_mm.txt；
#     排名库 *.genes_vs_motifs.rankings.feather（人 hg38、小鼠 mm10）；
#     motif 注释 motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl（小鼠用 .mgi-）
# 【输出结果】都在 results/08_scenic/
#   - tables/：08_adjacencies.tsv（阶段 I 的 TF-靶基因 关系）、08_motifs.csv、08_regulon_targets.csv（每个 regulon 的靶基因）、08_auc_matrix.csv（细胞 x regulon 活性，08c 的输入）
#   - objects/08_regulons.pkl
#   - logs/run_info_08b.json、pip_freeze_08b.txt
# 【怎么解读】
#   - regulon 名如 Gata1(+)：(+) = 激活型调控
#   - AUC 越高 = 该 TF 的靶基因在这个细胞里越集中高表达，即该 TF 越活跃
#   - regulon_targets 表可以查某个 TF 调控哪些基因，用于和差异基因、文献对照
# 【注意】
#   - 数据库版本要和 pyscenic 版本匹配：0.12.x 用 v10（clust）库，旧的 mc9nr 库对应 0.11.x
#   - 在 Spyder / Jupyter 里逐块运行（GRNBoost2 用 dask 多进程，直接 python 运行脚本在 macOS 上可能出错）；阶段 I 很耗时

# %%
# 步骤：导入包，设置路径和输出文件夹
# 输入：无
# 输出：路径变量；dirs（results/08_scenic/）
import glob
import os
import pickle
import pandas as pd
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
import fair_utils as fair

IN_FILE = "results/08_scenic/objects/08_expr_counts.csv"
DB_DIR = "data/cistarget"
TF_FILE = os.path.join(DB_DIR, "allTFs_hg38.txt")
DB_GLOB = os.path.join(DB_DIR, "*.genes_vs_motifs.rankings.feather")
MOTIF_FILE = os.path.join(DB_DIR, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
N_WORKERS = 4
SEED = 0
dirs = fair.setup_step("08_scenic")

# %%
# 步骤：读取表达矩阵、TF 列表和排名库
# 输入：IN_FILE、TF_FILE、DB_GLOB
# 输出：expr（细胞 x 基因 DataFrame）、tf_names、dbs
expr = pd.read_csv(IN_FILE, index_col=0)
tf_names = load_tf_names(TF_FILE)
dbs = [RankingDatabase(fname=f, name=os.path.basename(f)) for f in glob.glob(DB_GLOB)]
print(expr.shape, len(tf_names), dbs)

# %%
# 步骤：阶段 I — GRNBoost2 推断 TF → 靶基因 共表达关系
# 输入：expr、tf_names
# 输出：adjacencies（TF, target, importance）；tables/08_adjacencies.tsv
adjacencies = grnboost2(expression_data=expr, tf_names=tf_names, seed=SEED, verbose=True)
adjacencies.to_csv(fair.out(dirs, "tables", "adjacencies.tsv"), sep="\t", index=False)

# %%
# 步骤：阶段 II — 组装候选模块，用 motif 富集修剪成 regulon
# 输入：adjacencies、expr、dbs、MOTIF_FILE
# 输出：motifs、regulons；tables/08_motifs.csv、objects/08_regulons.pkl
modules = list(modules_from_adjacencies(adjacencies, expr))
motifs = prune2df(dbs, modules, MOTIF_FILE)
regulons = df2regulons(motifs)
motifs.to_csv(fair.out(dirs, "tables", "motifs.csv"))
with open(fair.out(dirs, "objects", "regulons.pkl"), "wb") as f:
    pickle.dump(regulons, f)

# %%
# 步骤：每个 regulon 的靶基因清单（便于查阅和引用）
# 输入：regulons
# 输出：tables/08_regulon_targets.csv（regulon, n_targets, targets）
regulon_targets = pd.DataFrame({
    "regulon": [r.name for r in regulons],
    "n_targets": [len(r.genes) for r in regulons],
    "targets": [";".join(sorted(r.genes)) for r in regulons],
})
regulon_targets.to_csv(fair.out(dirs, "tables", "regulon_targets.csv"), index=False)
regulon_targets.head()

# %%
# 步骤：阶段 III — AUCell 计算每个细胞的 regulon 活性
# 输入：expr、regulons
# 输出：auc_mtx（细胞 x regulon）；tables/08_auc_matrix.csv
auc_mtx = aucell(expr, regulons, num_workers=N_WORKERS)
auc_mtx.to_csv(fair.out(dirs, "tables", "auc_matrix.csv"))
auc_mtx.head()

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info_08b.json、logs/pip_freeze_08b.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()}, script="08b")
