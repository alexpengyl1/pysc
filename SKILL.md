---
name: pysc
description: Python/scanpy single-cell RNA-seq analysis (QC, integration, annotation, scCODA, GSEApy, pySCENIC, CellPhoneDB, LIANA, Palantir, scVelo, CellRank, inferCNVpy, Scirpy, scDRS, muon, metacells, CellOracle). Tells which conda environment and packages each step needs, the input data state it expects, and provides step-by-step templates.
---

# Python 单细胞转录组分析 / Python single-cell RNA-seq

模板来源：根据 Python 单细胞课程（ps2–ps30）提供的代码整理，已改成参数在顶部、按步骤拆分的脚本。
代码风格遵守 readable-coding skill：每块开头三行注释（步骤 / 输入 / 输出），一块只做一件事。
所有输出按 FAIR 原则组织（见「输出规范」）：每个步骤一个文件夹，图/表/对象/日志分开，文件名带步骤号，
自动记录参数和软件版本。

## 使用流程 / How to use this skill

用户要做某一步分析时，按顺序：

1. **确定步骤**：在下面的「步骤总表」找到对应模板和环境。
2. **先提示环境**：告诉用户要用哪个 conda 环境、怎么创建、需要哪些包（见「环境安装」），
   并让用户在该环境里运行 `python templates/00_check_env.py <环境名>`，确认没有 MISSING。
   用户的环境和这里不一样时，以用户导出的 `conda env export --no-builds` 为准。
3. **检查输入数据状态**：对照「数据约定」，确认上游文件存在、需要的是原始 counts 还是 log1p。
4. **复制模板和 `templates/fair_utils.py`** 到项目文件夹（两者放在一起，从项目根目录运行），
   只改顶部参数块（路径、物种、阈值、细胞类型名），说明选了哪些参数和原因。
   模板有两种格式，内容完全相同：`notebooks/*.ipynb`（Jupyter，用户默认用这个）和 `templates/*.py`（Spyder，`# %%` 分块）。
   notebook 已预设 kernel 名 = 环境名，打开后确认右上角 kernel 是对应环境。
5. 新增图或表时，按「输出规范」命名并保存，不要用 plt.show() 了事或把文件写到别处。
6. 结果解读看模板文件开头的「解读」段落；运行后每个步骤文件夹里的 README.md 列出了参数和全部输出。

## 输出规范（FAIR）/ Output convention

```
项目文件夹/
├── data/                       原始数据和参考文件（只读）
├── fair_utils.py               输出辅助模块（从 templates/ 复制）
├── 02_qc_cluster.ipynb …       分析 notebook（或 .py 脚本）
└── results/
    └── 02_qc_cluster/          每个步骤一个文件夹，名字 = 脚本名
        ├── README.md           自动生成：运行时间、conda 环境、参数、全部输出文件清单
        ├── figures/            02_qc_violin_by_sample.pdf + .png（PDF 矢量投稿用，PNG 300 dpi 预览用）
        ├── tables/             02_markers_leiden.csv
        ├── objects/            02_processed.h5ad
        └── logs/               run_info.json（时间、环境、平台、参数）、pip_freeze.txt（全部包版本）
```

- **文件名**：`<步骤号>_<内容>_<分组或基因>.<扩展名>`，只用小写英文、数字、下划线（不用空格和中文），
  例如 `05_umap_cell_type.pdf`、`07_deg_cd14_mono.csv`。文件被单独拷走也能看出来自哪一步、是什么。
- **代码写法**（fair_utils 提供，所有模板一致）：
  - 开头 `dirs = fair.setup_step("02_qc_cluster")` 建好四个子文件夹
  - 图：scanpy/scvelo 作图函数加 `show=False`，下一行 `fair.save_fig(dirs, "qc_violin_by_sample")`；
    返回图对象的（seaborn clustermap、plotnine）用 `fair.save_fig(dirs, "名字", fig=对象)`；
    CellRank 用自带 `save="14_xxx.pdf"`（`cr.settings.figdir` 已指向 figures/）
  - 表 / 对象路径：`fair.out(dirs, "tables", "markers_leiden.csv")` 自动加步骤号前缀
  - 结尾 `fair.write_run_info(dirs, params={k: v for k, v in globals().items() if k.isupper()})`
    收集顶部所有大写参数；一个步骤分几个脚本（08a/08b/08c、21a/21）时加 `script="08b"`，日志分开保存
- notebook 和脚本写法相同：在 Jupyter 里 `fair.save_fig` 会先存文件，再把图显示在 cell 下面。
- **参数全部写在脚本顶部、用大写变量名**，这样会被自动记录进 run_info.json 和 README.md。
- **格式**：表用 CSV，单细胞对象用 h5ad / h5mu（开放格式，R 和 Python 都能读），随机步骤固定 SEED。

## 步骤总表 / Step → template → environment

| 步骤 | 模板 | 环境 | 课程 |
|---|---|---|---|
| 读取、合并样本 | 01_read_data.py | sc-core | ps6 ps7 |
| QC、双细胞、标准化、HVG、PCA、聚类、marker | 02_qc_cluster.py | sc-core | ps2 |
| 整合 Harmony / BBKNN / ingest | 03_integrate_harmony_bbknn.py | sc-core | ps7 |
| 整合 scVI / scANVI + scib 评估 | 04_integrate_scvi.py | sc-scvi | ps10 |
| 注释：marker + CellTypist | 05_annotate.py | sc-core | ps2 ps9 |
| 细胞组成变化 scCODA | 06_composition_sccoda.py | sc-sccoda | ps11 |
| 差异基因 + GSEA / Enrichr | 07_de_gsea_gseapy.py | gseapy | ps12 |
| 转录调控 pySCENIC | 08a（sc-core）→ 08b（pyscenic）→ 08c（sc-core） | 见左 | ps13 |
| 细胞通讯 CellPhoneDB | 09_ccc_cellphonedb.py | cpdb | ps14 |
| 细胞通讯 LIANA+ | 10_ccc_liana.py | liana | ps15 |
| 拟时序 / 命运 Palantir | 11_trajectory_palantir.py | sc-palantir | ps16 |
| PHATE 可视化 | 12_embedding_phate.py | sc-core | ps17 |
| RNA velocity scVelo | 13_velocity_scvelo.py | sc-velo | ps18 |
| 命运概率 CellRank | 14_fate_cellrank.py | sc-palantir（拟时序）/ sc-velo（velocity） | ps19 |
| CNV inferCNVpy | 15_cnv_infercnvpy.py | sc-core | ps20 |
| 插补 MAGIC | 16_impute_magic.py | sc-core | ps21 |
| 药物靶点 drug2cell | 17_drug2cell.py | sc-core | ps22 |
| TCR/BCR Scirpy | 18_tcr_bcr_scirpy.py | sc-multi | ps23 |
| 基因模块 Hotspot | 19_modules_hotspot.py | sc-core | ps24 |
| GWAS 疾病相关细胞 scDRS | 20_disease_scdrs.py | sc-core | ps26 |
| 代谢通量 scFEA | 21a（sc-core 导出）→ 21_metabolism_scfea.sh（scFEA） | 见左 | ps27 |
| 多组学 RNA+ATAC muon | 22_multiome_muon.py | sc-multi | ps28 |
| 元细胞 metacells | 23_metacells.py | sc-metacell | ps29 |
| GRN CellOracle | 24_grn_celloracle.py | Docker | ps30 |
| 甲基化年龄 scAge（非转录组） | 25_methylation_scage.py | sc-core | ps25 |

上下游关系：01 → 02 → 03 或 04 → **05（05_annotated.h5ad 是所有下游的入口）** → 06–24 按需选择。
Palantir（11）或 scVelo（13）的结果再进入 CellRank（14）。

## 环境安装 / Environments

所有环境都可以用 `conda env create -f envs/<名字>.yml` 一步创建（yml 在本 skill 的 `envs/` 里）。
没有 yml 文件时用下面的命令。为什么分这么多环境：pySCENIC 需要 Python 3.7；scVI 依赖 PyTorch；
scCODA 依赖 TensorFlow；scvelo / liana / scFEA 课程里固定了旧版本。混装在一个环境里容易互相降级。

**sc-core**（主环境；对应课程的 spyder456）
```bash
conda create -n sc-core -c conda-forge python=3.10.14 scanpy=1.10.2 leidenalg python-igraph scikit-misc scikit-image pytables spyder-kernels=2.5 ipykernel
conda activate sc-core
pip install harmonypy==0.0.10 bbknn==1.6.0 celltypist phate scprep magic-impute infercnvpy drug2cell blitzgsea hotspotsc mplscience scdrs
```

**sc-scvi**
```bash
conda create -n sc-scvi -c conda-forge python=3.10.14 scanpy=1.10.2 leidenalg spyder-kernels=2.5 ipykernel
conda activate sc-scvi
pip install scvi-tools scib-metrics
```

**sc-sccoda**
```bash
conda create -n sc-sccoda -c conda-forge python=3.10.14 scanpy=1.10.2 spyder-kernels=2.5 ipykernel
conda activate sc-sccoda
pip install sccoda tf-keras pertpy
```

**gseapy**（Enrichr 基因集需要联网）
```bash
conda create -n gseapy -c conda-forge python=3.10.14 scanpy=1.10.2 lxml networkx spyder-kernels=2.5 spyder-notebook ipykernel
conda activate gseapy
pip install gseapy
```

**pyscenic**（Python 3.7）
```bash
conda create -n pyscenic -c conda-forge python=3.7 numpy cytoolz seaborn spyder-kernels ipykernel
conda activate pyscenic
pip install pyscenic
```
Apple Silicon（M 芯片）上 conda 可能没有原生 Python 3.7：在命令前加 `CONDA_SUBDIR=osx-64`（Rosetta 运行），
或改用官方 Docker 镜像 `aertslab/pyscenic`。cisTarget 数据库版本要和 pyscenic 版本匹配（0.12.x 用 v10 库）。

**cpdb**
```bash
conda create -n cpdb -c conda-forge python=3.10.14 spyder-kernels=2.5 spyder-notebook ipykernel
conda activate cpdb
pip install cellphonedb ktplotspy
```

**liana**
```bash
conda create -n liana -c conda-forge python=3.10.14 spyder-kernels=2.5 spyder-notebook ipykernel
conda activate liana
pip install "liana[common]" adjustText==0.8
```

**sc-palantir**（Palantir + CellRank 拟时序核）
```bash
conda create -n sc-palantir -c conda-forge python=3.10.14 scanpy=1.10.2 leidenalg spyder-kernels=2.5 ipykernel
conda activate sc-palantir
pip install palantir cellrank
# 可选：CellRank GAMR 基因趋势需要本机 R + mgcv，再 pip install rpy2
```

**sc-velo**（scVelo + CellRank velocity 核）
```bash
conda create -n sc-velo -c conda-forge python=3.10 scanpy=1.10.2 "numpy<2" leidenalg ipykernel notebook
conda activate sc-velo
pip install scvelo cellrank
# 生成 loom：pip install velocyto（需要 samtools），velocyto run10x -m repeat_msk.gtf <cellranger_dir> genes.gtf
```

**sc-multi**（muon + Scirpy）
```bash
conda create -n sc-multi -c conda-forge python=3.10.14 scanpy=1.10.2 leidenalg ipykernel notebook
conda activate sc-multi
pip install muon mudatasets pysam mofapy2 scirpy
```

**sc-metacell**
```bash
conda create -n sc-metacell -c conda-forge python=3.10.14 scanpy=1.10.2 seaborn ipykernel notebook
conda activate sc-metacell
pip install metacells
```

**scFEA**（命令行；先 `git clone https://github.com/changwn/scFEA`）
```bash
conda create -n scFEA -c pytorch -c conda-forge python=3.10 pytorch torchvision numpy pandas=1.5.3 scikit-learn matplotlib seaborn tqdm
conda activate scFEA
pip install magic-impute
python src/scFEA.py --help      # 在 scFEA 源码目录里检查
```

**CellOracle**（Docker，不用 conda）
```bash
docker pull kenjikamimoto126/celloracle_ubuntu:0.18.0
# 在项目根目录运行：把整个项目挂进容器，脚本、fair_utils.py 和 results/ 都能访问
docker run -dit --name celloracle_container -p 8889:8888 \
  -v "$PWD":/root/project -w /root/project \
  kenjikamimoto126/celloracle_ubuntu:0.18.0
docker exec -it celloracle_container bash
jupyter notebook --port=8888 --ip=0.0.0.0 --allow-root --no-browser   # 浏览器打开 http://localhost:8889
```
Apple Silicon 上这个镜像是 amd64 的，`docker run` 可能需要加 `--platform linux/amd64`。

**scAge**：用 sc-core 环境，`git clone https://github.com/alex-trapp/scAge`，脚本里 `sys.path.append` 指向该目录。

**在 Jupyter 里使用**：每个环境建好后注册成 kernel（名字和环境名一致，notebook 会自动匹配）：
```bash
conda activate sc-core
python -m ipykernel install --user --name sc-core --display-name sc-core
# 其他环境同理：sc-scvi、sc-sccoda、gseapy、pyscenic、cpdb、liana、sc-palantir、sc-velo、sc-multi、sc-metacell
```
Jupyter 本身只需装在一个地方（如 base 环境 `conda install -n base notebook`），在里面切换 kernel 即可。
CellOracle 的 notebook 在 Docker 容器的 Jupyter 里打开，提示找不到 kernel "celloracle" 时选容器里的 Python 3。

**在 Spyder 里使用**：该环境装好 `spyder-kernels`，Spyder → Preferences → Python interpreter 选这个环境的 python。

## 测试模板 / Testing

改了模板或新建/升级了环境后，在装 Jupyter 的环境里运行 `python tests/run_all_tests.py`（只测某几步：加步骤号，如 `05 07`），
Mac 上也可以双击 `tests/双击运行测试.command`（自动注册 kernel、装 nbclient、跑全部测试）。
它用 PBMC 3k 生成 4 个样本的测试数据（ctrl/stim 各 2 个，stim 里人为上调干扰素基因），按各 notebook 的 kernel 依次执行，
结果写入 `tests/test_report.md`（通过 / 失败 / 跳过 / 阻塞 + 报错）。需要外部大数据的步骤（cisTarget、loom、VDJ、multiome、GTF、
scAge 数据）缺文件时跳过，CellOracle 需在 Docker 里手动测。用户请你修 bug 时，先读 test_report.md 和 tests/test_project/executed/ 里出错的 notebook。

## 数据约定 / Data contract

**每个工具要求的矩阵状态不同，这是最容易出错的地方。**

步骤之间的交接文件：

| 文件 | 由谁生成 | 被谁读取 |
|---|---|---|
| data/sample_info.csv（sample, condition, …，用户准备） | 用户 | 01、06 |
| results/01_read_data/objects/01_merged_raw.h5ad | 01 | 02 |
| results/02_qc_cluster/objects/02_processed.h5ad | 02 | 03、04 |
| results/03_integrate/objects/03_integrated.h5ad | 03 | 05 |
| results/04_scvi/objects/04_scvi.h5ad | 04 | 05（替代 03 时改 IN_FILE） |
| **results/05_annotate/objects/05_annotated.h5ad** | 05 | 06–24 全部下游 |
| results/08_scenic/objects/08_expr_counts.csv | 08a | 08b |
| results/08_scenic/tables/08_auc_matrix.csv | 08b | 08c |
| results/11_palantir/objects/11_palantir.h5ad | 11 | 14 |
| results/13_scvelo/objects/13_velocity.h5ad | 13 | 14（velocity 核） |
| results/21_scfea/objects/21_expr.csv | 21a | 21 |
| results/24_celloracle/objects/24_input.h5ad | 24a（sc-core） | 24（Docker） |

05_annotated.h5ad（及 02、03）的固定结构：
- `X` = log1p（normalize_total 1e4 后取 log），**全部基因，未 scale**
- `layers["counts"]` = 原始 counts；`layers["log1p"]` = 同 X
- `obs`：`sample`、`data/sample_info.csv` 里的样本元数据（`condition` 等，01 读入时合并）、QC 指标、`leiden`、`cell_type`；`obsm`：`X_pca`、`X_umap`（整合后还有 `X_pca_harmony` / `X_scVI`）
- scale 只在 02 里对高变基因的临时副本做，用于 PCA，不保存

| 需要的输入 | 工具 |
|---|---|
| 原始 counts | scVI（layer="counts"）、HVG seurat_v3、Scrublet、Hotspot、scDRS、metacells、CellOracle（layers["raw_count"]）、pySCENIC（08a 导出 counts） |
| log1p、全基因、不 scale | CellTypist（要求每细胞总量 1e4）、差异分析、GSEApy、CellPhoneDB、LIANA（use_raw=False）、drug2cell（use_raw=False）、inferCNVpy、MAGIC、PHATE、scFEA（21a 导出） |
| spliced / unspliced | scVelo（loom，由 velocyto 生成） |
| 样本 x 细胞类型 计数 | scCODA（由 obs 汇总；需要生物学重复） |

## 物种 / Species

- 小鼠线粒体前缀 `mt-`，人 `MT-`（02、22、23 的参数）。
- **只支持人类基因名**：CellPhoneDB、drug2cell（ChEMBL）、scFEA 默认模块文件、VDJdb。小鼠数据要先用同源表
  （Ensembl BioMart 导出 mouse_symbol, human_symbol）转换，09 里有转换块。不要用简单的大写转换。
- 有小鼠选项：LIANA（`mouseconsensus`）、GSEApy（`organism="Mouse"`）、pySCENIC（mm10 库 + allTFs_mm）、
  CellOracle（`load_mouse_scATAC_atlas_base_GRN`）、scDRS（`load_gs` 同源转换）。
- CellTypist 模型大多是人类；小鼠肠道可用 `Adult_Mouse_Gut`。

## 已知问题 / Known pitfalls

- scanpy 1.10.2 + harmonypy ≥ 0.1：`harmony_integrate` 报 obsm shape 错误 → `harmonypy==0.0.10`。
- `sc.pp.scrublet` 需要 scikit-image。
- 重算 leiden 后再画 `rank_genes_groups_dotplot` 会报 dendrogram 类别数不一致 → 模板里用 `dendrogram=False`。
- BBKNN 新版报错时降到 1.6.0；liana 课程固定 adjustText 0.8；scFEA 需要 pandas 1.5.3。
- pySCENIC 的 GRNBoost2 用 dask 多进程：在 Spyder / Jupyter 里逐块运行。
- 课程代码里的演示性抽样（pySCENIC 只用 100 个细胞、Metacell 抽 5 万）不能用于正式分析；模板改为可设上限的参数。
- 一些包的作图函数内部会自己调用 plt.show()，在 Spyder/Jupyter 里 show 之后图就被关掉，
  `fair.save_fig` 会存成空白图。遇到空白图时：给该函数加 `show=False`，或改用它自带的 `save=` 参数存到 `dirs["figures"]`。
- 课程原代码的笔误已修正：`adjancencies`、`op_3_genes`、`scanpy.externalas`；
  过时 API 已替换：`AnnData.concatenate` → `ad.concat`，`rename_categories` → `map` 对照表。
- 逐细胞做条件间差异/富集会把细胞当作重复，p 值偏小；多个样本时考虑 pseudobulk，回复中提醒用户。

## 文件结构 / Files

```
pysc/
├── SKILL.md
├── envs/          # 12 个 conda yml（sc-core, sc-scvi, sc-sccoda, gseapy, pyscenic, cpdb, liana,
│                  #   sc-palantir, sc-velo, sc-multi, sc-metacell, scfea）
├── notebooks/     # 01–25 步骤模板的 Jupyter 版（.ipynb，kernel 预设为对应环境）
├── tests/         # make_test_data.py（PBMC 3k 测试数据）+ run_all_tests.py（按各自 kernel 跑全部 notebook，出 test_report.md）
└── templates/     # fair_utils.py（输出辅助）+ 00_check_env.py + 01–25 步骤模板（.py，Spyder 用）
```
.py 和 .ipynb 内容一致：notebook 由 .py 按 `# %%` 分块转换而来，文件开头说明变成 markdown，每个步骤前加小标题。
改了 .py 想同步到 notebook 时，可用 `jupytext --to ipynb 文件.py`（开头说明会留在代码注释里），或直接改 notebook。
如果当前环境里找不到 `templates/`，到用户的 GitHub 仓库 alexpengyl1/pysc 取。

## 附录：fair_utils.py / Appendix

找不到 `templates/fair_utils.py` 时，用下面的内容在项目文件夹里创建它：

```python
# FAIR 输出辅助模块：所有模板共用，放在和模板脚本同一个文件夹里
# 约定（每个步骤一个文件夹）：
#   results/<步骤>/figures/  图，PDF（矢量，投稿用）+ PNG（300 dpi，预览用）
#   results/<步骤>/tables/   表，CSV
#   results/<步骤>/objects/  数据对象，h5ad / h5mu / pickle
#   results/<步骤>/logs/     run_info.json（时间、环境、参数）、pip_freeze.txt（全部包版本）
#                            一个步骤分几个脚本时按脚本分开：run_info_08a.json、pip_freeze_08b.txt
#   results/<步骤>/README.md 自动生成：参数 + 输出文件清单
# 文件名：<步骤号>_<内容>_<分组或基因>.<扩展名>，只用小写英文、数字和下划线
# 兼容 Python 3.7（pyscenic 环境）

import datetime
import json
import os
import platform
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt

RESULTS_DIR = Path("results")
SUBDIRS = ("figures", "tables", "objects", "logs")


def setup_step(step):
    """建好 results/<step>/ 下的四个子文件夹，返回路径字典。"""
    step_dir = RESULTS_DIR / step
    dirs = {"step": step, "root": step_dir}
    for name in SUBDIRS:
        dirs[name] = step_dir / name
        dirs[name].mkdir(parents=True, exist_ok=True)
    return dirs


def out(dirs, kind, filename):
    """输出文件路径：results/<step>/<kind>/<步骤号>_<filename>。"""
    prefix = dirs["step"].split("_")[0]
    return dirs[kind] / f"{prefix}_{filename}"


def save_fig(dirs, name, fig=None):
    """保存当前图（或传入的图）为 PDF + PNG，然后显示并关闭。
    支持 matplotlib Figure、seaborn ClusterGrid、plotnine ggplot。"""
    if fig is None:
        fig = plt.gcf()
    pdf = out(dirs, "figures", f"{name}.pdf")
    png = out(dirs, "figures", f"{name}.png")
    if hasattr(fig, "savefig"):
        fig.savefig(pdf, bbox_inches="tight")
        fig.savefig(png, dpi=300, bbox_inches="tight")
        plt.show()
        plt.close("all")
    else:
        fig.save(str(pdf), verbose=False)
        fig.save(str(png), dpi=300, verbose=False)


def write_run_info(dirs, params, script=None):
    """记录运行信息和全部包版本，并生成该步骤的 README.md。
    一个步骤分几个脚本（如 08a/08b/08c）时传 script，日志按脚本分开保存。"""
    suffix = f"_{script}" if script else ""
    info = {
        "step": dirs["step"],
        "script": script or dirs["step"],
        "time": datetime.datetime.now().isoformat(timespec="seconds"),
        "conda_env": os.environ.get("CONDA_DEFAULT_ENV", ""),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "working_dir": str(Path.cwd()),
        "params": params,
    }
    with open(dirs["logs"] / f"run_info{suffix}.json", "w") as f:
        json.dump(info, f, indent=2, ensure_ascii=False, default=str)
    freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"], stdout=subprocess.PIPE, universal_newlines=True)
    (dirs["logs"] / f"pip_freeze{suffix}.txt").write_text(freeze.stdout)

    lines = [f"# {dirs['step']}", "",
             f"- 运行时间：{info['time']}",
             f"- conda 环境：{info['conda_env']}（Python {info['python']}）",
             f"- 完整包版本：logs/pip_freeze{suffix}.txt", "", f"## 参数（{info['script']}）", ""]
    lines += [f"- `{key}` = {value}" for key, value in params.items()]
    for kind in ("figures", "tables", "objects", "logs"):
        files = sorted(p.name for p in dirs[kind].iterdir() if p.is_file())
        lines += ["", f"## {kind}", ""] + [f"- {kind}/{name}" for name in files]
    (dirs["root"] / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("已写入", dirs["root"] / "README.md")
```
