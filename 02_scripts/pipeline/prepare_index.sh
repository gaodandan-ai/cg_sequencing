#!/bin/bash
#SBATCH --job-name=prepare_index
#SBATCH --partition=qcpu_23a
#SBATCH --error=08_logs/index_build/prepare_index_%j.err
#SBATCH --output=08_logs/index_build/prepare_index_%j.out
#SBATCH -n 8  # 优化：索引构建最多用8核，足够且不浪费
#SBATCH --mem=8G  # 优化：谷棒基因组小，8G内存完全够用
# 如需邮件提醒，可在本地提交脚本中自行添加 --mail-type / --mail-user

# ============================================================================
# 谷氨酸棒杆菌参考基因组索引准备脚本
# 功能：构建 BWA 和 Samtools 索引
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# 可通过环境变量覆盖工具路径，例如：BWA=/path/to/bwa bash 02_scripts/pipeline/prepare_index.sh
bwa="${BWA:-bwa}"
samtools="${SAMTOOLS:-samtools}"

# 参考基因组和索引目录路径
reference="${REFERENCE:-"${PROJECT_ROOT}/00_reference/genome.fna"}"
index_dir="${INDEX_DIR:-"${PROJECT_ROOT}/00_reference/index"}"
# Slurm日志目录（提前创建，避免提交失败）
log_dir="${LOG_DIR:-"${PROJECT_ROOT}/08_logs/index_build"}"

# 打印任务信息
echo "==========================================================================="
echo "                    开始执行谷氨酸棒杆菌参考基因组索引构建"
echo "==========================================================================="
echo "参考基因组: $reference"
echo "索引输出目录: $index_dir"
echo "构建内容: BWA索引 + Samtools索引"
echo "==========================================================================="

# 关键：提前创建索引目录和日志目录（确保存在）
mkdir -p "$index_dir"
mkdir -p "$log_dir"

# Step 1: 构建/验证BWA索引
echo -e "\nStep 1: 验证/构建BWA索引"
echo "--------------------------------------------------------------------------"
# 检查index_dir下的BWA索引是否存在（因为指定了-p ${index_dir}/genome）
bwa_index_check="${index_dir}/genome.bwt"
if [ -f "$bwa_index_check" ] && [ -f "${index_dir}/genome.sa" ]; then
    echo "✅  BWA索引已存在，无需重复构建"
else
    echo "正在构建BWA索引..."
    "$bwa" index -p "${index_dir}/genome" "$reference"
    echo "✅  BWA索引构建成功（输出到${index_dir}）"
fi

# Step 2: 构建/验证Samtools索引（补全else分支和fi闭合，修复语法错误）
echo -e "\nStep 2: 验证/构建Samtools索引"
echo "--------------------------------------------------------------------------"
samtools_index_check="${reference}.fai"
if [ -f "$samtools_index_check" ]; then
    echo "✅  Samtools索引已存在，无需重复构建"
else
    echo "正在构建Samtools索引..."
    "$samtools" faidx "$reference"
    echo "✅  Samtools索引构建成功（输出到${reference}.fai）"
fi

# Step 3: 最终验证
echo -e "\n==========================================================================="
echo "                    索引构建流程执行完成！"
echo "---------------------------------------------------------------------------"
echo "✅  BWA索引位置: ${index_dir}/genome.*"
echo "✅  Samtools索引位置: ${reference}.fai"
echo "==========================================================================="

# 脚本正常退出
exit 0
