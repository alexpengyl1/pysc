#!/bin/zsh -i
# 双击即可运行全部模板测试。结果：tests/test_report.md、tests/test_console.log
cd "$(dirname "$0")"

# 步骤：找到 conda（终端里已能用 conda 时跳过）
if ! command -v conda >/dev/null 2>&1; then
  for d in ~/miniconda3 ~/anaconda3 ~/miniforge3 ~/mambaforge /opt/miniconda3 /opt/anaconda3 /opt/homebrew/Caskroom/miniforge/base; do
    [ -f "$d/etc/profile.d/conda.sh" ] && source "$d/etc/profile.d/conda.sh" && break
  done
fi
command -v conda >/dev/null 2>&1 || { echo "找不到 conda，请告诉 Claude 你的 conda 装在哪里"; read; exit 1; }

# 步骤：把已建好的环境注册成同名 kernel（没建的环境会列出来，对应步骤会报“缺 kernel”）
for env in sc-core sc-scvi sc-sccoda gseapy pyscenic cpdb liana sc-palantir sc-velo sc-multi sc-metacell; do
  if conda env list | grep -q "^$env "; then
    conda run -n "$env" python -m ipykernel install --user --name "$env" --display-name "$env" >/dev/null 2>&1 \
      && echo "kernel 已注册：$env" || echo "注册失败（环境里可能没装 ipykernel）：$env"
  else
    echo "环境不存在：$env"
  fi
done

# 步骤：在 base 环境里装测试工具并运行全部测试
conda activate base
python -m pip install -q nbclient nbformat ipykernel
python run_all_tests.py 2>&1 | tee test_console.log

echo ""
echo "测试结束，可以关掉这个窗口，回到 Claude 说一声“跑完了”。"
read
