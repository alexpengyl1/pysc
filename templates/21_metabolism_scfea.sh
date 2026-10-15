# 21 scFEA: infer metabolic flux in each cell from gene expression
#
# [Purpose]
#   scFEA uses a graph neural network to estimate, from the expression of metabolic genes, the flux of each metabolic module (e.g. glycolysis, TCA cycle) in each cell,
#   and the accumulation or consumption of each metabolite. This shell script calls the scFEA source code.
# [Pipeline position] Upstream 21a_export_for_scfea.py -> this step (end point)
# [Environment] Run 21a in sc-core first; then, in the scFEA environment, run this script from the project root:
#   conda activate scFEA && bash 21_metabolism_scfea.sh (course ps27)
#   Setup: git clone https://github.com/changwn/scFEA and put the path in SCFEA_DIR
# [Input]
#   - results/21_scfea/objects/21_expr.csv (genes x cells, log-normalised, exported by 21a)
#   - scFEA's module gene table and stoichiometry matrix (human: module_gene_m168.csv + cmMat_c70_m168.csv; mouse: the corresponding files with mouse in scFEA/data/)
# [Output] All in results/21_scfea/
#   - tables/: flux matrix (cells x metabolic modules) and balance matrix (cells x metabolites), csv
#   - logs/: run_info_21.txt (time, parameters, command), pip_freeze_21.txt
# [How to read the results]
#   - The flux matrix can be used as a new feature matrix: with tables/21_cell_info.csv, compare glycolysis, TCA and other module fluxes between cell types
#   - The balance matrix reflects accumulation (positive) or consumption (negative) pressure on each metabolite
# [Notes]
#   - scFEA trains a model on every run, so many cells make it slow

set -e
SCFEA_DIR="$HOME/software/scFEA"
STEP_DIR="$PWD/results/21_scfea"
MODULE_FILE="module_gene_m168.csv"
STOICHIOMETRY_FILE="cmMat_c70_m168.csv"
IMPUTATION="True"

# Step: record run time, parameters and package versions
# Input: the variables above
# Output: logs/run_info_21.txt, logs/pip_freeze_21.txt
mkdir -p "$STEP_DIR/tables" "$STEP_DIR/logs"
echo "time: $(date '+%Y-%m-%dT%H:%M:%S')" > "$STEP_DIR/logs/run_info_21.txt"
echo "conda_env: $CONDA_DEFAULT_ENV" >> "$STEP_DIR/logs/run_info_21.txt"
echo "MODULE_FILE=$MODULE_FILE STOICHIOMETRY_FILE=$STOICHIOMETRY_FILE IMPUTATION=$IMPUTATION" >> "$STEP_DIR/logs/run_info_21.txt"
python -m pip freeze > "$STEP_DIR/logs/pip_freeze_21.txt"

# Step: run scFEA (imputation recommended for 10x data)
# Input: objects/21_expr.csv, module gene table, stoichiometry matrix
# Output: flux and balance matrices (csv) in tables/
cd "$SCFEA_DIR"
python src/scFEA.py \
    --data_dir data \
    --input_dir "$STEP_DIR/objects" \
    --test_file 21_expr.csv \
    --moduleGene_file "$MODULE_FILE" \
    --stoichiometry_matrix "$STOICHIOMETRY_FILE" \
    --res_dir "$STEP_DIR/tables" \
    --sc_imputation "$IMPUTATION"

# Step: list the output files
# Input: tables/
# Output: printed file list
ls -lh "$STEP_DIR/tables"
