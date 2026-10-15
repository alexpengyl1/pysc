# FAIR output helpers shared by all templates; keep this file in the same folder as the template scripts
# Convention (one folder per step):
#   results/<step>/figures/  figures, PDF (vector, for papers) + PNG (300 dpi, for preview)
#   results/<step>/tables/   tables, CSV
#   results/<step>/objects/  data objects, h5ad / h5mu / pickle
#   results/<step>/logs/     run_info.json (time, environment, parameters), pip_freeze.txt (all package versions)
#                            a step split into several scripts keeps separate logs: run_info_08a.json, pip_freeze_08b.txt
#   results/<step>/README.md written automatically: parameters + list of output files
# File names: <step number>_<content>_<group or gene>.<extension>, lowercase letters, digits and underscores only
# Compatible with Python 3.7 (pyscenic environment)

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
    """Create the four subfolders under results/<step>/ and return a dict of paths."""
    step_dir = RESULTS_DIR / step
    dirs = {"step": step, "root": step_dir}
    for name in SUBDIRS:
        dirs[name] = step_dir / name
        dirs[name].mkdir(parents=True, exist_ok=True)
    return dirs


def out(dirs, kind, filename):
    """Output path: results/<step>/<kind>/<step number>_<filename>."""
    prefix = dirs["step"].split("_")[0]
    return dirs[kind] / f"{prefix}_{filename}"


def save_fig(dirs, name, fig=None):
    """Save the current figure (or the one passed in) as PDF + PNG, then show and close it.
    Supports matplotlib Figure, seaborn ClusterGrid and plotnine ggplot."""
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
    """Record run information and all package versions, and write this step's README.md.
    When a step is split into several scripts (e.g. 08a/08b/08c), pass script so logs are kept separately."""
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
             f"- Run time: {info['time']}",
             f"- Conda environment: {info['conda_env']} (Python {info['python']})",
             f"- Full package versions: logs/pip_freeze{suffix}.txt", "", f"## Parameters ({info['script']})", ""]
    lines += [f"- `{key}` = {value}" for key, value in params.items()]
    for kind in ("figures", "tables", "objects", "logs"):
        files = sorted(p.name for p in dirs[kind].iterdir() if p.is_file())
        lines += ["", f"## {kind}", ""] + [f"- {kind}/{name}" for name in files]
    (dirs["root"] / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("Written:", dirs["root"] / "README.md")
