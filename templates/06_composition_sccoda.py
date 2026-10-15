# 06 Cell-type composition changes (scCODA)
#
# [Purpose]
#   Tests whether the proportion of each cell type really changes between conditions (e.g. irradiated vs control).
#   Proportions are compositional (one type going up pushes the others down), so scCODA uses a Bayesian model with a reference cell type.
#   Workflow: sample x cell type counts -> stacked bar plot and box plots -> choose reference cell type -> fit model and sample -> credible changes
# [Pipeline position] Upstream 05 -> this step (end point)
# [Environment] sc-sccoda (TensorFlow; course ps11)
# [Input]
#   - results/05_annotate/objects/05_annotated.h5ad (obs['cell_type'], obs['sample'])
#   - data/sample_info.csv (sample, condition; one row per sample)
# [Output] All in results/06_sccoda/
#   - tables/: 06_counts_sample_celltype.csv (counts table), 06_effects.csv, 06_intercepts.csv, 06_credible_effects.csv
#   - figures/: 06_stacked_barplot_condition, 06_boxplot_condition, 06_reference_selection
#   - objects/06_sccoda_result.pkl: model result, can be reloaded
# [How to read the results]
#   - Stacked bar / box plots: a visual view of proportions per group; each point in the box plot is one sample
#   - effects table: Final Parameter not equal to 0 = credible change; log2-fold change is relative to the control group
#   - credible_effects: use this for conclusions; an apparent "decrease" in other cell types may just be compositional squeezing
#   - reference_selection plot: choose a reference present in every sample with a stable proportion (high presence, low dispersion); the result depends on this choice
# [Notes]
#   - Each condition needs at least 2 biological replicate samples; the unit of analysis is the sample, not the cell
#   - The default FDR of 0.05 is strict; 0.2 is more sensitive but gives more false positives

# %%
# Step: import packages, set parameters and output folder
# Input: none
# Output: parameter variables; dirs (results/06_sccoda/)
import pandas as pd
import scanpy as sc
from sccoda.util import cell_composition_data as dat
from sccoda.util import comp_ana as mod
from sccoda.util import data_visualization as viz
import fair_utils as fair

IN_FILE = "results/05_annotate/objects/05_annotated.h5ad"
SAMPLE_INFO = "data/sample_info.csv"
CELLTYPE_KEY = "cell_type"
SAMPLE_KEY = "sample"
CONDITION_KEY = "condition"
REFERENCE = "automatic"          # or a cell type name
FDR = 0.05
dirs = fair.setup_step("06_sccoda")

# %%
# Step: read the annotated data and the sample table
# Input: IN_FILE, SAMPLE_INFO
# Output: adata, sample_info (index = sample name)
adata = sc.read_h5ad(IN_FILE)
sample_info = pd.read_csv(SAMPLE_INFO, index_col="sample")
sample_info

# %%
# Step: build and save the sample x cell type count table
# Input: adata.obs, sample_info
# Output: counts (AnnData: rows = samples, columns = cell types); tables/06_counts_sample_celltype.csv
counts = dat.from_scanpy(adata, cell_type_identifier=CELLTYPE_KEY, sample_identifier=SAMPLE_KEY,
                         covariate_df=sample_info)
counts.to_df().to_csv(fair.out(dirs, "tables", "counts_sample_celltype.csv"))
counts

# %%
# Step: stacked bar plot of cell composition per condition
# Input: counts
# Output: figures/06_stacked_barplot_condition.pdf/.png
viz.stacked_barplot(counts, feature_name=CONDITION_KEY)
fair.save_fig(dirs, "stacked_barplot_condition")

# %%
# Step: box plots of each cell type's proportion by condition
# Input: counts
# Output: figures/06_boxplot_condition.pdf/.png
viz.boxplots(counts, feature_name=CONDITION_KEY, add_dots=True)
fair.save_fig(dirs, "boxplot_condition")

# %%
# Step: reference cell type selection plot (high presence, low dispersion suits a reference)
# Input: counts
# Output: figures/06_reference_selection.pdf/.png
viz.rel_abundance_dispersion_plot(data=counts, abundant_threshold=0.9)
fair.save_fig(dirs, "reference_selection")

# %%
# Step: build the model and run HMC sampling
# Input: counts
# Output: result (CAResult)
model = mod.CompositionalAnalysis(counts, formula=CONDITION_KEY, reference_cell_type=REFERENCE)
result = model.sample_hmc()
result.set_fdr(est_fdr=FDR)
result.summary()

# %%
# Step: save effects, intercepts and credible effects tables
# Input: result
# Output: tables/06_effects.csv, 06_intercepts.csv, 06_credible_effects.csv
result.effect_df.to_csv(fair.out(dirs, "tables", "effects.csv"))
result.intercept_df.to_csv(fair.out(dirs, "tables", "intercepts.csv"))
result.credible_effects().to_csv(fair.out(dirs, "tables", "credible_effects.csv"))

# %%
# Step: save the result object
# Input: result
# Output: objects/06_sccoda_result.pkl
result.save(str(fair.out(dirs, "objects", "sccoda_result.pkl")))

# %%
# Step: record parameters and package versions, write this step's README
# Input: parameters at the top of this script (UPPERCASE variables)
# Output: logs/run_info.json, logs/pip_freeze.txt, README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
