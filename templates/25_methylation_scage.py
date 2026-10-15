# 25 scAge：单细胞 DNA 甲基化表观遗传年龄（不是转录组分析）
#
# 【这份代码做什么】
#   用 bulk 甲基化数据为每个 CpG 位点建立“甲基化水平 ~ 年龄”的线性模型作为参考，
#   再根据单细胞甲基化测序的结果，预测每个细胞的表观遗传年龄。
#   流程：读取 bulk 训练矩阵 → 每个 CpG 拟合年龄模型 → 单细胞覆盖文件二值化 → 预测年龄 → 按分组作图
# 【在流程中的位置】独立分析（单细胞甲基化数据，不依赖 01–05）
# 【运行环境】sc-core（需要 pytables）；scAge 不是 pip 包，要 git clone https://github.com/alex-trapp/scAge（课程 ps25）
# 【输入数据】
#   - bulk 甲基化训练矩阵 .h5：样本 x CpG，含 Age 列（TRAIN_MATRIX）
#   - 单细胞 Bismark 覆盖文件 *.cov.gz（放在 COV_DIR）
# 【输出结果】都在 results/25_scage/
#   - objects/：25_reference.tsv（每个 CpG 的 甲基化~年龄 模型）、25_sc_processed/（二值化后的单细胞甲基化）
#   - tables/：scAge 预测结果（每个细胞的预测年龄）、25_training_samples.csv
#   - figures/25_predicted_age_by_group
# 【怎么解读】
#   - 每个细胞的预测表观年龄；箱线图比较不同组（如年轻 vs 老年、处理 vs 对照）的年龄差异
#   - 组织特异参考模型更准，多组织模型适合来源不明的细胞
#   - CPG_PERCENTILE = 1 表示只用与年龄相关性最强的前 1% 位点（单组织模型）；多组织模型可用 10

# %%
# 步骤：导入包，设置路径和输出文件夹
# 输入：无
# 输出：参数变量；dirs（results/25_scage/）
import glob
import multiprocessing
import sys
import pandas as pd
import matplotlib.pyplot as plt
import fair_utils as fair

SCAGE_DIR = "software/scAge"        # 放 scAge.py 的目录
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
# 步骤：读取 bulk 训练矩阵，保存训练样本信息
# 输入：TRAIN_MATRIX
# 输出：train（样本 x CpG + Age 等列）；tables/25_training_samples.csv
train = pd.read_hdf(TRAIN_MATRIX)
train[[c for c in train.columns if not c.startswith("chr")]].to_csv(fair.out(dirs, "tables", "training_samples.csv"))
train.shape

# %%
# 步骤：为每个 CpG 拟合 甲基化 ~ 年龄 线性模型，作为参考
# 输入：train
# 输出：objects/25_reference.tsv
reference_file = str(fair.out(dirs, "objects", "reference.tsv"))
scAge.construct_reference(training_DNAm_matrix=train, output_path=reference_file, n_cores=N_CORES, chunksize=200)

# %%
# 步骤：单细胞覆盖文件二值化（甲基化 >= 50% 记 1）
# 输入：COV_DIR 下的 *.cov.gz
# 输出：objects/25_sc_processed/
processed_dir = str(fair.out(dirs, "objects", "sc_processed")) + "/"
scAge.process_coverage(cov_directory=COV_DIR, output_path=processed_dir,
                       n_cores=N_CORES, max_met=100, split="_", chunksize=1, binarization="round")

# %%
# 步骤：预测每个细胞的表观年龄
# 输入：processed_dir、reference_file
# 输出：tables/ 下的预测结果文件
scAge.run_scAge(single_cell_dir_or_dict=processed_dir, single_cell_set_name=SET_NAME,
                reference_data=reference_file, output_path=str(dirs["tables"]) + "/",
                selection_mode="percentile", CpG_parameter=CPG_PERCENTILE,
                min_age=MIN_AGE, max_age=MAX_AGE, age_step=AGE_STEP, n_cores=N_CORES)

# %%
# 步骤：读取预测结果，按细胞分组画预测年龄（分组规则按自己的细胞命名改；列名以 scAge 实际输出为准）
# 输入：tables/ 下的预测结果
# 输出：predictions；figures/25_predicted_age_by_group.pdf/.png
predictions = pd.read_csv(glob.glob(str(dirs["tables"] / f"*{SET_NAME}*"))[0], sep="\t", index_col=0)
predictions["group"] = predictions.index.str.split("_").str[0]
predictions.boxplot(column="PredictedAge", by="group")
plt.ylabel("Predicted epigenetic age (months)")
fair.save_fig(dirs, "predicted_age_by_group")

# %%
# 步骤：记录参数、软件版本，生成本步骤 README
# 输入：本脚本顶部的参数（全大写变量）
# 输出：logs/run_info.json、logs/pip_freeze.txt、README.md
fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})
