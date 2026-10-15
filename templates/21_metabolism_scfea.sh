# 21 scFEA：从基因表达推断每个细胞的代谢通量
#
# 【这份代码做什么】
#   scFEA 用图神经网络，根据代谢相关基因的表达估计每个细胞中各代谢模块（如糖酵解、TCA 循环）的通量，
#   以及各代谢物的积累或消耗。这是一个命令行脚本，调用 scFEA 源码运行。
# 【在流程中的位置】上游 21a_export_for_scfea.py → 本步（终点分析）
# 【运行环境】先在 sc-core 里运行 21a；再在 scFEA 环境里，从项目根目录运行本脚本：
#   conda activate scFEA && bash 21_metabolism_scfea.sh（课程 ps27）
#   准备：git clone https://github.com/changwn/scFEA，把路径填到 SCFEA_DIR
# 【输入数据】
#   - results/21_scfea/objects/21_expr.csv（基因 x 细胞，log 标准化表达，21a 导出）
#   - scFEA 自带的模块基因表和化学计量矩阵（人：module_gene_m168.csv + cmMat_c70_m168.csv；小鼠用 scFEA/data/ 里带 mouse 的对应文件）
# 【输出结果】都在 results/21_scfea/
#   - tables/：通量矩阵（细胞 x 代谢模块）和平衡矩阵（细胞 x 代谢物），csv
#   - logs/：run_info_21.txt（时间、参数、命令）、pip_freeze_21.txt
# 【怎么解读】
#   - 通量矩阵可以当成新的特征矩阵：结合 tables/21_cell_info.csv 按细胞类型比较糖酵解、TCA 等模块的通量
#   - 平衡矩阵反映代谢物的积累（正）或消耗（负）压力
# 【注意】
#   - 每次运行都要训练模型，细胞多时很慢

set -e
SCFEA_DIR="$HOME/software/scFEA"
STEP_DIR="$PWD/results/21_scfea"
MODULE_FILE="module_gene_m168.csv"
STOICHIOMETRY_FILE="cmMat_c70_m168.csv"
IMPUTATION="True"

# 步骤：记录运行时间、参数和包版本
# 输入：上面的变量
# 输出：logs/run_info_21.txt、logs/pip_freeze_21.txt
mkdir -p "$STEP_DIR/tables" "$STEP_DIR/logs"
echo "time: $(date '+%Y-%m-%dT%H:%M:%S')" > "$STEP_DIR/logs/run_info_21.txt"
echo "conda_env: $CONDA_DEFAULT_ENV" >> "$STEP_DIR/logs/run_info_21.txt"
echo "MODULE_FILE=$MODULE_FILE STOICHIOMETRY_FILE=$STOICHIOMETRY_FILE IMPUTATION=$IMPUTATION" >> "$STEP_DIR/logs/run_info_21.txt"
python -m pip freeze > "$STEP_DIR/logs/pip_freeze_21.txt"

# 步骤：运行 scFEA（10x 数据建议插补）
# 输入：objects/21_expr.csv、模块基因表、化学计量矩阵
# 输出：tables/ 下的通量和平衡矩阵 csv
cd "$SCFEA_DIR"
python src/scFEA.py \
    --data_dir data \
    --input_dir "$STEP_DIR/objects" \
    --test_file 21_expr.csv \
    --moduleGene_file "$MODULE_FILE" \
    --stoichiometry_matrix "$STOICHIOMETRY_FILE" \
    --res_dir "$STEP_DIR/tables" \
    --sc_imputation "$IMPUTATION"

# 步骤：查看输出文件
# 输入：tables/
# 输出：打印的文件列表
ls -lh "$STEP_DIR/tables"
