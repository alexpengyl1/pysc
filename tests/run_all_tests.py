# Test every notebook template in one go: run each from scratch on small test data in its own conda kernel and write a test report
#
# Usage (Mac terminal, in the environment where Jupyter is installed, usually base):
#   conda activate base
#   pip install nbclient nbformat          # first time only
#   cd "<pysc folder>/tests"
#   python run_all_tests.py                # all steps
#   python run_all_tests.py 05 07 09       # selected steps only (their upstream results must already exist, e.g. retesting after a bug fix)
#
# Requirement: each environment registered as a Jupyter kernel of the same name (see "Using Jupyter" in SKILL.md)
# Output (all in tests/):
#   test_project/            test project: data/ (test data), results/ (FAIR outputs of each step), executed/ (executed notebooks)
#   test_report.md           pass / fail / skip / blocked for each step, run time, error and failing step for failures
#   test_report.json         same, machine-readable
# Needs internet: downloads PBMC 3k, CellTypist models, Enrichr gene sets, the CellPhoneDB database

import datetime
import json
import re
import shutil
import sys
import time
import traceback
from pathlib import Path

import nbformat
from jupyter_client.kernelspec import KernelSpecManager
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

TESTS_DIR = Path(__file__).resolve().parent
SKILL_DIR = TESTS_DIR.parent
PROJECT = TESTS_DIR / "test_project"
CELL_TIMEOUT = 3600          # maximum run time of a single cell in seconds

# Test order, upstream dependencies and skip conditions (skipped when an external file is missing)
STEPS = [
    ("make_test_data", [], None),
    ("01_read_data", ["make_test_data"], None),
    ("02_qc_cluster", ["01_read_data"], None),
    ("03_integrate_harmony_bbknn", ["02_qc_cluster"], None),
    ("04_integrate_scvi", ["02_qc_cluster"], None),
    ("05_annotate", ["03_integrate_harmony_bbknn"], None),
    ("06_composition_sccoda", ["05_annotate"], None),
    ("07_de_gsea_gseapy", ["05_annotate"], None),
    ("08a_export_for_scenic", ["05_annotate"], None),
    ("08b_pyscenic", ["08a_export_for_scenic"], "data/cistarget/*.feather"),
    ("08c_scenic_plot", ["08b_pyscenic"], "data/cistarget/*.feather"),
    ("09_ccc_cellphonedb", ["05_annotate"], None),
    ("10_ccc_liana", ["05_annotate"], None),
    ("11_trajectory_palantir", ["05_annotate"], None),
    ("12_embedding_phate", ["05_annotate"], None),
    ("13_velocity_scvelo", ["05_annotate"], "data/velocyto/*.loom"),
    ("14_fate_cellrank", ["11_trajectory_palantir"], None),
    ("15_cnv_infercnvpy", ["05_annotate"], "data/genes.gtf"),
    ("16_impute_magic", ["05_annotate"], None),
    ("17_drug2cell", ["05_annotate"], None),
    ("18_tcr_bcr_scirpy", ["05_annotate"], "data/S1/filtered_contig_annotations.csv"),
    ("19_modules_hotspot", ["05_annotate"], None),
    ("20_disease_scdrs", ["05_annotate"], None),
    ("21a_export_for_scfea", ["05_annotate"], None),
    ("22_multiome_muon", [], "data/multiome/filtered_feature_bc_matrix.h5"),
    ("23_metacells", ["05_annotate"], None),
    ("24_grn_celloracle", ["05_annotate"], "DOCKER"),
    ("25_methylation_scage", [], "data/scage/sc_data_raw/*"),
]

# Parameters replaced for testing (a whole line or a whole multi-line dict/list) and code snippets
FIRST_TYPE = 'sc.read_h5ad(IN_FILE).obs.query("cell_type == \'{}\'").index[0]'
PARAMS = {
    "01_read_data": {"SAMPLES": '{s: RAW_DIR / s for s in ["S1", "S2", "S3", "S4"]}'},
    "05_annotate": {"MARKERS": '{"T": ["CD3D", "IL7R"], "CD8_T": ["CD8A"], "B": ["MS4A1", "CD79A"], '
                               '"NK": ["GNLY", "NKG7"], "Mono": ["CD14", "LYZ", "FCGR3A"], "DC": ["FCER1A"], "Mk": ["PPBP"]}'},
    "07_de_gsea_gseapy": {"CELLTYPE": '"CD14_Mono"', "CASE, CONTROL": '"stim", "ctrl"'},
    "09_ccc_cellphonedb": {"ITERATIONS": "100", "SOURCE": '"CD14_Mono"', "TARGETS": '"CD4_T|B"',
                           "GENES": '["CD74", "MIF", "CCL5", "HLA-DRA"]'},
    "10_ccc_liana": {"N_PERMS": "100", "SOURCES": '["CD14_Mono", "B"]', "TARGETS": '["CD4_T", "NK"]'},
    "11_trajectory_palantir": {
        "LINEAGE_TYPES": '["CD14_Mono", "FCGR3A_Mono", "DC"]',
        "START_CELL": FIRST_TYPE.format("CD14_Mono"),
        "TERMINAL_STATES": "{" + FIRST_TYPE.format("FCGR3A_Mono") + ': "FCGR3A_Mono", ' + FIRST_TYPE.format("DC") + ': "DC"}',
        "N_WAYPOINTS": "200",
        "GENES": '["CD14", "LYZ", "FCGR3A", "FCER1A"]'},
    "14_fate_cellrank": {"START_CLUSTER": '"CD14_Mono"', "N_STATES": "3", "GENES": '["CD14", "LYZ", "FCGR3A", "FCER1A"]'},
    "15_cnv_infercnvpy": {"REFERENCE_TYPES": '["CD4_T", "CD8_T", "NK"]', "TUMOR_CLUSTERS": '["0"]'},
    "16_impute_magic": {"GENES": '["CD14", "LYZ", "FCGR3A"]'},
    "19_modules_hotspot": {"CELLTYPE": '"CD4_T"', "N_GENES": "200"},
    "20_disease_scdrs": {"N_CTRL": "200"},
}
CODE = {
    "04_integrate_scvi": [("model.train()", "model.train(max_epochs=20)")],
    "05_annotate": [('adata.obs["cell_type"] = adata.obs[CLUSTER_KEY].map(annotation).astype("category")',
                     'test_types = pd.read_csv("data/test_cell_types.csv", index_col="barcode")["cell_type"]  # for testing\n'
                     'adata.obs["cell_type"] = test_types.reindex(adata.obs_names.str.rsplit("-", n=1).str[0]).values\n'
                     'adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")')],
}


def set_param(source, name, value):
    """Replace `NAME = ...` (including a multi-line dict/list) with `NAME = value`."""
    lines = source.split("\n")
    for i, line in enumerate(lines):
        if re.match(rf"^{re.escape(name)}\s*=", line):
            end = i
            if line.split("#")[0].rstrip().endswith(("{", "[", "(")):
                while not lines[end].strip().startswith(("}", "]", ")")):
                    end += 1
            lines[i:end + 1] = [f"{name} = {value}   # replaced for testing"]
            return "\n".join(lines), True
    return source, False


def make_test_data_notebook():
    """Wrap make_test_data.py as a notebook and run it in the sc-core kernel."""
    source = (TESTS_DIR / "make_test_data.py").read_text(encoding="utf-8")
    cells = [nbformat.v4.new_code_cell(c.strip()) for c in re.split(r"^# %%[^\n]*\n", source, flags=re.M) if c.strip()]
    nb = nbformat.v4.new_notebook(cells=cells)
    nb.metadata["kernelspec"] = {"name": "sc-core", "display_name": "sc-core", "language": "python"}
    return nb


def prepare(step):
    """Read a notebook, apply test parameters and code replacements, return (notebook, replacements not found)."""
    if step == "make_test_data":
        return make_test_data_notebook(), []
    nb = nbformat.read(SKILL_DIR / "notebooks" / f"{step}.ipynb", as_version=4)
    missing = []
    for name, value in PARAMS.get(step, {}).items():
        done = False
        for cell in nb.cells:
            if cell.cell_type == "code" and not done:
                cell.source, done = set_param(cell.source, name, value)
        if not done:
            missing.append(f"parameter {name}")
    for old, new in CODE.get(step, []):
        hits = [c for c in nb.cells if c.cell_type == "code" and old in c.source]
        for cell in hits:
            cell.source = cell.source.replace(old, new)
        if not hits:
            missing.append(f"code {old[:40]}")
    return nb, missing


def failing_step_title(nb):
    """Find the '# Step:' description of the cell that failed."""
    for cell in nb.cells:
        if cell.cell_type == "code" and any(o.get("output_type") == "error" for o in cell.get("outputs", [])):
            match = re.search(r"^# Step: (.+)$", cell.source, flags=re.M)
            return match.group(1) if match else cell.source.splitlines()[0]
    return ""


def run_step(step, kernels):
    nb, missing = prepare(step)
    kernel = nb.metadata["kernelspec"]["name"]
    if kernel not in kernels:
        return {"status": "FAIL", "kernel": kernel, "seconds": 0,
                "error": f"No Jupyter kernel named {kernel}. Run conda activate {kernel}, then "
                         f"python -m ipykernel install --user --name {kernel} --display-name {kernel}"}
    start = time.time()
    client = NotebookClient(nb, timeout=CELL_TIMEOUT, kernel_name=kernel, resources={"metadata": {"path": str(PROJECT)}})
    result = {"kernel": kernel, "not_applied": missing}
    try:
        client.execute()
        result["status"] = "PASS"
    except CellExecutionError as error:
        message = re.sub(r"\x1b\[[0-9;]*m", "", str(error))      # strip terminal colour codes
        result.update(status="FAIL", failed_at=failing_step_title(nb),
                      error="\n".join(message.strip().splitlines()[-12:]))
    except Exception:
        result.update(status="FAIL", error=traceback.format_exc()[-1500:])
    result["seconds"] = round(time.time() - start)
    nbformat.write(nb, PROJECT / "executed" / f"{step}.ipynb")
    return result


def write_report(results):
    icons = {"PASS": "✅ pass", "FAIL": "❌ fail", "SKIP": "⏭ skip", "BLOCKED": "⛔ blocked"}
    lines = ["# Template test report", "", f"- Time: {datetime.datetime.now():%Y-%m-%d %H:%M}",
             f"- Test project: {PROJECT}", "", "| Step | kernel | Result | Time (s) | Notes |", "|---|---|---|---|---|"]
    for step, r in results.items():
        note = r.get("failed_at") or r.get("reason") or ""
        lines.append(f"| {step} | {r.get('kernel', '')} | {icons[r['status']]} | {r.get('seconds', '')} | {note} |")
    for step, r in results.items():
        if r["status"] == "FAIL" or r.get("not_applied"):
            lines += ["", f"## {step}", ""]
            if r.get("not_applied"):
                lines.append(f"Test replacements not found: {r['not_applied']}")
            if r.get("error"):
                lines += ["```", r["error"], "```"]
    (TESTS_DIR / "test_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (TESTS_DIR / "test_report.json").write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")


# Main: prepare the test project folder, run steps in order, update the report after each step
(PROJECT / "executed").mkdir(parents=True, exist_ok=True)
shutil.copy(SKILL_DIR / "templates" / "fair_utils.py", PROJECT / "fair_utils.py")
kernels = KernelSpecManager().find_kernel_specs()
selected = sys.argv[1:]
report_file = TESTS_DIR / "test_report.json"
results = json.loads(report_file.read_text(encoding="utf-8")) if selected and report_file.exists() else {}

for step, deps, needs in STEPS:
    if selected and not any(step.startswith(s) for s in selected):
        continue
    blocked = [d for d in deps if results.get(d, {}).get("status") in ("FAIL", "BLOCKED", "SKIP")]
    if needs == "DOCKER":
        results[step] = {"status": "SKIP", "reason": "run by hand inside the CellOracle Docker container"}
    elif needs and not list(PROJECT.glob(needs)):
        results[step] = {"status": "SKIP", "reason": f"missing external data {needs}"}
    elif blocked:
        results[step] = {"status": "BLOCKED", "reason": f"upstream {', '.join(blocked)} did not pass"}
    else:
        print(f"Running {step} ...", flush=True)
        results[step] = run_step(step, kernels)
    print(f"  {step}: {results[step]['status']}", flush=True)
    write_report(results)

print("Report:", TESTS_DIR / "test_report.md")
