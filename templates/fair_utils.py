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
