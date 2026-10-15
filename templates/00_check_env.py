# 00 Environment check
#
# [Purpose]
#   Checks, one by one, whether the packages an environment needs are installed, and prints their versions.
#   Run it after creating or upgrading an environment so you do not discover a missing package halfway through an analysis.
# [Environment] The environment being checked, e.g. conda activate sc-core, then python 00_check_env.py sc-core
# [Input]
#   - Command-line argument: environment name (sc-core, sc-scvi, sc-sccoda, gseapy, pyscenic, cpdb, liana, sc-palantir, sc-velo, sc-multi, sc-metacell, scFEA)
# [Output]
#   - Printed to screen: Python version + version of each package; no files are written
# [How to read the results]
#   - A version number = installed; MISSING = not installed, install it following "Environments" in SKILL.md or envs/<env>.yml
#   - If a version differs from the one pinned in SKILL.md (e.g. harmonypy, bbknn, adjustText), later steps may fail; switch to the pinned version first

import sys
from importlib.metadata import version, PackageNotFoundError

# Step: packages to check for each environment
# Input: none
# Output: PACKAGES (environment name -> list of packages)
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

# Step: print the Python version and each package version; missing packages show MISSING
# Input: command-line argument = environment name
# Output: printed list of versions
env_name = sys.argv[1]
print("python", sys.version.split()[0])
for package in PACKAGES[env_name]:
    try:
        print(package, version(package))
    except PackageNotFoundError:
        print(package, "MISSING")
