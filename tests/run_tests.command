#!/bin/zsh -i
# Double-click to run all template tests. Results: tests/test_report.md, tests/test_console.log
cd "$(dirname "$0")"

# Step: find conda (skipped when conda already works in the terminal)
if ! command -v conda >/dev/null 2>&1; then
  for d in ~/miniconda3 ~/anaconda3 ~/miniforge3 ~/mambaforge /opt/miniconda3 /opt/anaconda3 /opt/homebrew/Caskroom/miniforge/base; do
    [ -f "$d/etc/profile.d/conda.sh" ] && source "$d/etc/profile.d/conda.sh" && break
  done
fi
command -v conda >/dev/null 2>&1 || { echo "conda not found; tell Claude where conda is installed"; read; exit 1; }

# Step: register existing environments as Jupyter kernels of the same name (missing environments are listed; their steps will report a missing kernel)
for env in sc-core sc-scvi sc-sccoda gseapy pyscenic cpdb liana sc-palantir sc-velo sc-multi sc-metacell; do
  if conda env list | grep -q "^$env "; then
    conda run -n "$env" python -m ipykernel install --user --name "$env" --display-name "$env" >/dev/null 2>&1 \
      && echo "kernel registered: $env" || echo "registration failed (ipykernel may be missing): $env"
  else
    echo "environment not found: $env"
  fi
done

# Step: install the test tools in base and run all tests
conda activate base
python -m pip install -q nbclient nbformat ipykernel
python run_all_tests.py 2>&1 | tee test_console.log

echo ""
echo "Tests finished. You can close this window and tell Claude they are done."
read
