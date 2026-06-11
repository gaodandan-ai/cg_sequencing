#!/bin/bash
# ============================================================================
# 谷氨酸棒杆菌重测序分析 - 一键流水线
# 依次执行：索引构建 → 全流程分析
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# 颜色输出
green() { echo -e "\033[32m$1\033[0m"; }
red()   { echo -e "\033[31m$1\033[0m"; }
cyan()  { echo -e "\033[36m$1\033[0m"; }

echo "============================================================================"
cyan "                   谷氨酸棒杆菌重测序分析流水线"
echo "============================================================================"
echo "项目目录: ${PROJECT_ROOT}"
echo "开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================================"
echo ""

# --------------------------------------------------
# Step 1: 检查输入
# --------------------------------------------------
if [ ! -f "${PROJECT_ROOT}/00_reference/genome.fna" ]; then
    red "❌ 错误：参考基因组不存在 → ${PROJECT_ROOT}/00_reference/genome.fna"
    echo "请将参考基因组 FASTA 文件放入 00_reference/，并命名为 genome.fna"
    exit 1
fi

fastq_count=$(find "${PROJECT_ROOT}/01_raw_data" -maxdepth 1 \( -name "*.fq.gz" -o -name "*.fastq.gz" \) 2>/dev/null | wc -l)
if [ "$fastq_count" -eq 0 ]; then
    red "❌ 错误：01_raw_data/ 下没有 FASTQ 文件"
    echo "请将双端测序 FASTQ 文件放入 01_raw_data/"
    exit 1
fi
echo "✅ 输入检查通过：参考基因组 + ${fastq_count} 个 FASTQ 文件"
echo ""

# --------------------------------------------------
# Step 2: 检查 Conda 环境
# --------------------------------------------------
if [ -z "${CONDA_DEFAULT_ENV:-}" ]; then
    echo "⚠️  未检测到 Conda 环境激活"
    echo "   建议先运行: conda activate cg_seq"
    echo "   继续执行，但可能因工具缺失而失败..."
    echo ""
fi

# --------------------------------------------------
# Step 3: 构建参考基因组索引
# --------------------------------------------------
green "═══ Step 1/2: 构建参考基因组索引 ═══"
echo ""

# 如果索引已存在则跳过
if [ -f "${PROJECT_ROOT}/00_reference/index/genome.bwt" ] && \
   [ -f "${PROJECT_ROOT}/00_reference/genome.fna.fai" ]; then
    echo "✅ 索引已存在，跳过索引构建步骤"
else
    bash "${SCRIPT_DIR}/CG_PrepareIndex.sh"
    echo ""
fi

echo ""

# --------------------------------------------------
# Step 4: 运行主流程
# --------------------------------------------------
green "═══ Step 2/2: 运行重测序分析主流程 ═══"
echo ""

mkdir -p "${PROJECT_ROOT}/08_logs/sequencing_analysis"
log_file="${PROJECT_ROOT}/08_logs/sequencing_analysis/pipeline_$(date '+%Y%m%d_%H%M%S').log"

bash "${SCRIPT_DIR}/CG_GenomeReseq.sh" 2>&1 | tee "$log_file"
pipeline_exit=${PIPESTATUS[0]}

echo ""

# --------------------------------------------------
# Step 5: 完成
# --------------------------------------------------
echo "============================================================================"
if [ "$pipeline_exit" -eq 0 ]; then
    green "✅ 全流程执行成功！"
else
    red "❌ 流水线执行失败（退出码: $pipeline_exit）"
    echo "   详细日志: $log_file"
fi
echo "结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "日志文件: $log_file"
echo "============================================================================"

exit "$pipeline_exit"
