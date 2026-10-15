# 25 scAge: epigenetic age from single-cell DNA methylation (not RNA-seq)
#
# [Purpose]
#   Fits a "methylation level ~ age" linear model for each CpG from bulk methylation data as a reference,
#   then predicts the epigenetic age of each cell from single-cell methylation sequencing.
#   Workflow: read bulk training matrix -> fit age model per CpG -> binarise single-cell coverage files -> predict age -> plot by group
# [Pipeline position] Standalone (single-cell methylation data, does not depend on 01-05)
# [Environment] sc-core (needs pytables); scAge is not a pip package: git clone https://github.com/alex-trapp/scAge (course ps25)
# [Input]
#   - Bulk methylation training matrix .h5: samples x CpGs, with an Age column (TRAIN_MATRIX)
#   - Single-cell Bismark coverage files *.cov.gz (in COV_DIR)
# [Output] All in results/25_scage/
#   - objects/: 25_reference.tsv (methylation ~ age model per CpG), 25_sc_processed/ (binarised single-cell methylation)
#   - tables/: scAge predictions (predicted age per cell), 25_training_samples.csv
#   - figures/25_predicted_age_by_group
# [How to read the results]
#   - Predicted epigenetic age per cell; the box plot compares groups (e.g. young vs old, treated vs control)
#   - Tissue-specific reference models are more accurate; multi-tissue models suit cells of unknown origin
#   - CPG_PERCENTILE = 1 uses only the top 1% of sites most correlated with age (single-tissue model); multi-tissue models can use 10

# %%
# Step: import packages, set paths and output folder
# Input: none
# Output: parameter variables; dirs (results/25_scage/)
import glob
import multiprocessing
import sys
import pandas as pd
import matplotlib.pyplot as plt
import fair_utils as fair

SCAGE_DIR = "software/scAge"        # folder containing scAge.py
sys.path.append(SCAGE_DIR)
import scAge

TRAIN_MATRIX = "data/scage/Thompson_Liver_BL6_bulk_matrix.h5"
COV_DIR = "data/scage/sc_data_raw/"
SET_NAME = "my_cells"
CPG_PERCENTILE = 1
MIN_AGE, MAX_AGE, AGE_STEP = -20, 60, 0.1
N_CORES = multiprocessing.cpu_count()
dirs = fair.setup_step("25_scage")

# %%
# Step: read the bulk training matrix, save the training sample information
# Input: TRAIN_MATRIX
# Output: train (samples x CpGs + Age and other columns); tables/25_training_samples.csv
train = pd.read_hdf(TRAIN_MATRIX)
train[[c for c in train.columns if not c.startswith("chr")]].to_csv(fair.out(dirs, "tables", "training_samples.csv"))
train.shape

# %%
# Step: fit a methylation ~ age linear model for each CpG as the reference
# Input: train
# Output: objects/25_reference.tsv
reference_file = str(fair.out(dirs, "objects", "reference.tsv"))
scAge.construct_reference(training_DNAm_matrix=train, output_path=reference_file, n_cores=N_CORES, chunksize=200)

# %%
# Step: binarise single-cell coverage files (methylation >= 50% counts as 1)
# Input: *.cov.gz in COV_DIR
# Output: objects/25_sc_processed/
processed_dir = str(fair.out(dirs, "objects", "sc_processed")) + "/"
scAge.process_coverage(cov_directory=COV_DIR, output_path=processed_dir,
                       n_cores=N_CORES, max_met=100, split="_", chunksize=1, binarization="round")

# %%
# Step: predict the epigenetic age of each cell
# Input: processed_dir, reference_file
# Output: prediction files in tables/
scAge.run_scAge(single_cell_dir_or_dict=processed_dir, single_cell_set_name=SET_NAME,
                reference_data=reference_file, output_path=str(dirs["tables"]) + "/",
                selection_mode="percentile", CpG_parameter=CPG_PERCENTILE,
                min_age=MIN_AGE, max_age=MAX_AGE, age_step=AGE_STEP, n_cores=N_CORES)

# %%
# Step: read the predictions and plot predicted age by group (adapt the grouping to your cell naming; column names follow the actual scAge output)
# Input: predictions in tables/
# Output: predictions; figures/25_predicted_age_by_group.pdf/.png
predictions = pd.read_csv(glob.glob(str(dirs["tables"] / f"*{SET_NAME}*"))[0], sep="\t", index_col=0)
predictions["group"] = predictions.index.str.split("_").str[0]
predictions.boxplot(column="PredictedAge", by="group")
plt.ylabel("Predicted epigenetic age (months)")
fair.save_fig(dirs, "predicted_age_by_group")

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
