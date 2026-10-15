# 08b pySCENIC: infer transcription factor regulons and their activity in each cell
#
# [Purpose]
#   Finds which target genes each transcription factor (TF) regulates and scores the activity of these regulons in every cell. Three stages:
#     Stage I   GRNBoost2: infers TF -> target links from co-expression
#     Stage II  cisTarget: keeps only links where the TF motif is enriched near the target promoter, giving regulons
#     Stage III AUCell: scores each regulon's activity (AUC) in each cell
# [Pipeline position] Upstream 08a -> this step -> downstream 08c_scenic_plot
# [Environment] pyscenic (Python 3.7; course ps13)
# [Input]
#   - results/08_scenic/objects/08_expr_counts.csv (exported by 08a)
#   - Databases in data/cistarget/ (https://resources.aertslab.org/cistarget/):
#     TF list allTFs_hg38.txt / allTFs_mm.txt;
#     ranking databases *.genes_vs_motifs.rankings.feather (human hg38, mouse mm10);
#     motif annotation motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl (mouse: .mgi-)
# [Output] All in results/08_scenic/
#   - tables/: 08_adjacencies.tsv (stage I TF-target links), 08_motifs.csv, 08_regulon_targets.csv (targets of each regulon), 08_auc_matrix.csv (cells x regulon activity, input to 08c)
#   - objects/08_regulons.pkl
#   - logs/run_info_08b.json, pip_freeze_08b.txt
# [How to read the results]
#   - Regulon names such as Gata1(+): (+) = activating
#   - Higher AUC = the TF's targets are more concentrated among the highly expressed genes of that cell, i.e. the TF is more active
#   - The regulon_targets table lists the genes each TF regulates, for comparison with DE genes and the literature
# [Notes]
#   - Database versions must match pyscenic: 0.12.x uses the v10 (clust) databases, the older mc9nr databases go with 0.11.x
#   - Run block by block in Spyder / Jupyter (GRNBoost2 uses dask multiprocessing; running the script directly with python can fail on macOS); stage I is slow

# %%
# Step: import packages, set paths and output folder
# Input: none
# Output: path variables; dirs (results/08_scenic/)
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
# Step: read the expression matrix, TF list and ranking databases
# Input: IN_FILE, TF_FILE, DB_GLOB
# Output: expr (cells x genes DataFrame), tf_names, dbs
expr = pd.read_csv(IN_FILE, index_col=0)
tf_names = load_tf_names(TF_FILE)
dbs = [RankingDatabase(fname=f, name=os.path.basename(f)) for f in glob.glob(DB_GLOB)]
print(expr.shape, len(tf_names), dbs)

# %%
# Step: stage I - GRNBoost2 infers TF -> target co-expression links
# Input: expr, tf_names
# Output: adjacencies (TF, target, importance); tables/08_adjacencies.tsv
adjacencies = grnboost2(expression_data=expr, tf_names=tf_names, seed=SEED, verbose=True)
adjacencies.to_csv(fair.out(dirs, "tables", "adjacencies.tsv"), sep="\t", index=False)

# %%
# Step: stage II - build candidate modules and prune them into regulons by motif enrichment
# Input: adjacencies, expr, dbs, MOTIF_FILE
# Output: motifs, regulons; tables/08_motifs.csv, objects/08_regulons.pkl
modules = list(modules_from_adjacencies(adjacencies, expr))
motifs = prune2df(dbs, modules, MOTIF_FILE)
regulons = df2regulons(motifs)
motifs.to_csv(fair.out(dirs, "tables", "motifs.csv"))
with open(fair.out(dirs, "objects", "regulons.pkl"), "wb") as f:
    pickle.dump(regulons, f)

# %%
# Step: target gene list of each regulon (for lookup and citation)
# Input: regulons
# Output: tables/08_regulon_targets.csv (regulon, n_targets, targets)
regulon_targets = pd.DataFrame({
    "regulon": [r.name for r in regulons],
    "n_targets": [len(r.genes) for r in regulons],
    "targets": [";".join(sorted(r.genes)) for r in regulons],
})
regulon_targets.to_csv(fair.out(dirs, "tables", "regulon_targets.csv"), index=False)
regulon_targets.head()

# %%
# Step: stage III - AUCell scores regulon activity in each cell
# Input: expr, regulons
# Output: auc_mtx (cells x regulons); tables/08_auc_matrix.csv
auc_mtx = aucell(expr, regulons, num_workers=N_WORKERS)
auc_mtx.to_csv(fair.out(dirs, "tables", "auc_matrix.csv"))
auc_mtx.head()

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info_08b.json, logs/pip_freeze_08b.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()}, script="08b")
