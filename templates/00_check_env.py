# 00 环境检查
#
# 【这份代码做什么】
#   在某个 conda 环境里逐个检查该环境需要的包是否已安装，并打印版本号。
#   每次新建或升级环境后先跑一次，避免分析跑到一半才发现缺包。
# 【运行环境】要检查的那个环境本身，例如 conda activate sc-core 后运行 python 00_check_env.py sc-core
# 【输入数据】
#   - 命令行参数：环境名（sc-core、sc-scvi、sc-sccoda、gseapy、pyscenic、cpdb、liana、sc-palantir、sc-velo、sc-multi、sc-metacell、scFEA）
# 【输出结果】
#   - 屏幕打印：Python 版本 + 每个包的版本，不生成文件
# 【怎么解读】
#   - 显示版本号 = 已装好；显示 MISSING = 没装，按 SKILL.md「环境安装」或 envs/<环境>.yml 补装
#   - 版本和 SKILL.md 里固定的版本不同时（如 harmonypy、bbknn、adjustText），后续步骤可能报错，优先改成固定版本

import sys
from importlib.metadata import version, PackageNotFoundError

# 步骤：每个环境需要检查的包
# 输入：无
# 输出：PACKAGES（环境名 -> 包名列表）
PACKAGES = {
    "sc-core": ["scanpy", "anndata", "leidenalg", "scikit-misc", "harmonypy", "bbknn", "celltypist",
                "phate", "scprep", "magic-impute", "infercnvpy", "drug2cell", "blitzgsea",
                "hotspotsc", "scdrs", "tables"],
    "sc-scvi": ["scanpy", "scvi-tools", "torch", "scib-metrics"],
    "sc-sccoda": ["scanpy", "sccoda", "tensorflow", "tf-keras", "pertpy"],
    "gseapy": ["scanpy", "gseapy", "networkx"],
    "pyscenic": ["pyscenic", "arboreto", "ctxcore", "pandas", "numpy"],
    "cpdb": ["cellphonedb", "anndata", "ktplotspy"],
    "liana": ["liana", "scanpy", "adjustText"],
    "sc-palantir": ["scanpy", "palantir", "cellrank"],
    "sc-velo": ["scanpy", "scvelo", "cellrank", "numpy"],
    "sc-multi": ["scanpy", "muon", "mofapy2", "pysam", "scirpy"],
    "sc-metacell": ["scanpy", "metacells"],
    "scFEA": ["torch", "pandas", "magic-impute"],
}

# 步骤：打印 Python 版本和每个包的版本，没装的显示 MISSING
# 输入：命令行参数 = 环境名
# 输出：打印的版本列表
env_name = sys.argv[1]
print("python", sys.version.split()[0])
for package in PACKAGES[env_name]:
    try:
        print(package, version(package))
    except PackageNotFoundError:
        print(package, "MISSING")
