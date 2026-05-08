#!/bin/bash
# ============================================================================
# 谷氨酸棒杆菌重测序分析脚本（双输出版，升级稳健版）
# 支持原始双端命名：
#   1) xxx_1.fastq.gz / xxx_2.fastq.gz
#   2) xxx_R1.fastq.gz / xxx_R2.fastq.gz
#   3) xxx_R1_001.fastq.gz / xxx_R2_001.fastq.gz
#   4) xxx_1.fq.gz / xxx_2.fq.gz
#   5) xxx_R1.fq.gz / xxx_R2.fq.gz
#   6) xxx_R1_001.fq.gz / xxx_R2_001.fq.gz
#
# 流程：
# 原始质控 → 修剪 → 修剪后质控 → 比对 → 排序 → 去重 → 变异检测
# → 基础质量过滤（不按AF硬过滤）
# → 按AF拆分为 polymorphic / highfreq
# → 两套结果分别注释与统计
# ============================================================================

set -euo pipefail
shopt -s nullglob

# -------------------------- 1. 环境与工具配置 --------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

if [ -f "${PROJECT_ROOT}/config.env" ]; then
    # shellcheck disable=SC1091
    source "${PROJECT_ROOT}/config.env"
fi

export JAVA_OPTS="${JAVA_OPTS:-"-Xmx16g -Xms4g"}"

bwa="${BWA:-bwa}"
trim_galore="${TRIM_GALORE:-trim_galore}"
fastqc="${FASTQC:-fastqc}"
multiqc="${MULTIQC:-multiqc}"
samtools="${SAMTOOLS:-samtools}"
bcftools="${BCFTOOLS:-bcftools}"
cutadapt="${CUTADAPT:-cutadapt}"
snpeff="${SNPEFF:-snpEff}"
picard="${PICARD:-picard}"
python3_bin="${PYTHON3:-python3}"

# -------------------------- 2. 参数设置 --------------------------
reference="${REFERENCE:-"${PROJECT_ROOT}/00_reference/genome.fna"}"
bwa_index="${BWA_INDEX:-"${PROJECT_ROOT}/00_reference/index/genome"}"

SNPEFF_DATA_DIR="${SNPEFF_DATA_DIR:-"${PROJECT_ROOT}/07_annotation"}"
CG_DB_NAME="${CG_DB_NAME:-Corynebacterium_glutamicum_ATCC13032}"
if [ -z "${SNPEFF_CONFIG:-}" ]; then
    if [ -n "${CONDA_PREFIX:-}" ]; then
        snpeff_configs=( "${CONDA_PREFIX}"/share/snpeff*/snpEff.config )
        SNPEFF_CONFIG="${snpeff_configs[0]:-}"
    else
        SNPEFF_CONFIG=""
    fi
fi

fastq_dir="${FASTQ_DIR:-"${PROJECT_ROOT}/01_raw_data"}"
fastqc_raw="${FASTQC_RAW:-"${PROJECT_ROOT}/03_qc/raw"}"
fastqc_clean="${FASTQC_CLEAN:-"${PROJECT_ROOT}/03_qc/trimmed"}"
trimgalore_dir="${CLEAN_DATA_DIR:-"${PROJECT_ROOT}/04_clean_data"}"
alignment="${ALIGNMENT_DIR:-"${PROJECT_ROOT}/05_alignment/sam"}"
sorted="${SORTED_BAM_DIR:-"${PROJECT_ROOT}/05_alignment/sorted_bam"}"
dedup="${DEDUP_BAM_DIR:-"${PROJECT_ROOT}/05_alignment/dedup_bam"}"

variant_raw="${VARIANT_RAW_DIR:-"${PROJECT_ROOT}/06_variant_calling/raw"}"
variant_qc="${VARIANT_QC_DIR:-"${PROJECT_ROOT}/06_variant_calling/quality_filtered"}"
variant_poly="${VARIANT_POLY_DIR:-"${PROJECT_ROOT}/06_variant_calling/polymorphic"}"
variant_high="${VARIANT_HIGH_DIR:-"${PROJECT_ROOT}/06_variant_calling/highfreq"}"

variant_poly_annot="${VARIANT_POLY_ANNOT_DIR:-"${PROJECT_ROOT}/06_variant_calling/annotated/polymorphic"}"
variant_high_annot="${VARIANT_HIGH_ANNOT_DIR:-"${PROJECT_ROOT}/06_variant_calling/annotated/highfreq"}"

variant_tables="${VARIANT_TABLES_DIR:-"${PROJECT_ROOT}/06_variant_calling/tables"}"
helper_dir="${HELPER_DIR:-"${PROJECT_ROOT}/02_scripts/helpers"}"
log_dir="${LOG_DIR:-"${PROJECT_ROOT}/08_logs/sequencing_analysis"}"

make_abs_path() {
    case "$1" in
        /*) printf '%s\n' "$1" ;;
        *) printf '%s/%s\n' "$PROJECT_ROOT" "$1" ;;
    esac
}

reference="$(make_abs_path "$reference")"
bwa_index="$(make_abs_path "$bwa_index")"

for dir_var in fastq_dir fastqc_raw fastqc_clean trimgalore_dir alignment sorted dedup \
    variant_raw variant_qc variant_poly variant_high variant_poly_annot variant_high_annot \
    variant_tables helper_dir log_dir SNPEFF_DATA_DIR; do
    value="${!dir_var}"
    value="$(make_abs_path "$value")"
    value="${value%/}/"
    printf -v "$dir_var" '%s' "$value"
done

# 阈值参数
THREADS="${THREADS:-4}"
SORT_MEM="${SORT_MEM:-1G}"
MIN_DP_CALL="${MIN_DP_CALL:-10}"
MIN_DP_AF="${MIN_DP_AF:-20}"
MIN_AF_POLY="${MIN_AF_POLY:-0.01}"
MAX_AF_POLY="${MAX_AF_POLY:-0.99}"
MIN_AF_HIGH="${MIN_AF_HIGH:-0.80}"
MIN_QUAL="${MIN_QUAL:-30}"

# -------------------------- 3. 初始化 --------------------------
echo "=== 初始化：创建目录 & 检查依赖 ==="

mkdir -p \
  "$fastqc_raw" "$fastqc_clean" "$trimgalore_dir" \
  "$alignment" "$sorted" "$dedup" \
  "$variant_raw" "$variant_qc" "$variant_poly" "$variant_high" \
  "$variant_poly_annot" "$variant_high_annot" \
  "$variant_tables" "$helper_dir" "$log_dir"

tools=("$bwa" "$trim_galore" "$fastqc" "$multiqc" "$samtools" "$bcftools" "$cutadapt" "$snpeff" "$picard" "$python3_bin")
for tool in "${tools[@]}"; do
    if ! command -v "$tool" >/dev/null 2>&1 && [ ! -x "$tool" ]; then
        echo "❌ 错误：工具不存在或不可执行 → $tool"
        exit 1
    fi
done

if [ ! -f "$reference" ]; then
    echo "❌ 错误：参考基因组不存在 → $reference"
    exit 1
fi

if [ ! -f "${bwa_index}.bwt" ]; then
    echo "⚠️ BWA索引不存在，自动构建..."
    "$bwa" index -p "$bwa_index" "$reference"
fi

if [ ! -d "${SNPEFF_DATA_DIR}/${CG_DB_NAME}" ] || [ ! -f "${SNPEFF_DATA_DIR}/${CG_DB_NAME}/snpEffectPredictor.bin" ]; then
    echo "❌ 错误：snpEff自定义数据库不存在！"
    echo "请先构建数据库："
    echo "$snpeff build -Xmx16g -gff3 -v -noCheckCds -noCheckProtein -dataDir $SNPEFF_DATA_DIR $CG_DB_NAME"
    exit 1
else
    echo "✅ snpEff数据库检查通过：$CG_DB_NAME"
fi

snpeff_config_args=()
if [ -n "$SNPEFF_CONFIG" ] && [ -f "$SNPEFF_CONFIG" ]; then
    snpeff_config_args=( -config "$SNPEFF_CONFIG" )
elif [ -n "$SNPEFF_CONFIG" ]; then
    echo "⚠️ 未找到 SNPEFF_CONFIG 指定文件：$SNPEFF_CONFIG，将使用 snpEff 默认配置"
fi

if find "$fastq_dir" -maxdepth 1 \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | grep -q .; then
    fastq_n=$(find "$fastq_dir" -maxdepth 1 \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | wc -l)
    echo "✅ 原始数据检查通过：找到 ${fastq_n} 个 FASTQ 文件"
else
    echo "❌ 错误：$fastq_dir 下没有 .fastq.gz 或 .fq.gz 文件"
    exit 1
fi

# -------------------------- 4. 检查辅助Python脚本 --------------------------

split_py="${helper_dir}/split_vcf_by_af.py"
pair_py="${helper_dir}/detect_fastq_pairs.py"

for helper in "$split_py" "$pair_py"; do
    if [ ! -f "$helper" ]; then
        echo "❌ 错误：辅助脚本不存在 → $helper"
        exit 1
    fi
    chmod +x "$helper"
done

# -------------------------- 5. 原始数据FastQC --------------------------
echo -e "\n=== Step 1: 原始数据FastQC质控 ==="
mapfile -t raw_fastqs < <(find "$fastq_dir" -maxdepth 1 \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | sort)

if [ "${#raw_fastqs[@]}" -eq 0 ]; then
    echo "❌ 错误：未找到原始FASTQ文件"
    exit 1
fi

"$fastqc" -t "$THREADS" -o "$fastqc_raw" "${raw_fastqs[@]}"
"$multiqc" -o "$fastqc_raw" "$fastqc_raw"
echo "✅ 原始数据质控完成"

# -------------------------- 6. 识别双端配对 --------------------------
echo -e "\n=== Step 2: 识别原始FASTQ双端配对 ==="
pair_tsv="${trimgalore_dir}/fastq_pairs.tsv"
"$python3_bin" "$pair_py" --indir "$fastq_dir" --out "$pair_tsv"

sample_count=$(tail -n +2 "$pair_tsv" | wc -l)
echo "✅ 已识别到 $sample_count 个双端样本"

# -------------------------- 7. Trim Galore --------------------------
echo -e "\n=== Step 3: Trim Galore! 修剪 ==="

tail -n +2 "$pair_tsv" | while IFS=$'\t' read -r sample fq1 fq2; do
    trimmed_done=false
    for r1_candidate in \
        "${trimgalore_dir}/${sample}_1_val_1.fq.gz" \
        "${trimgalore_dir}/${sample}_R1_val_1.fq.gz" \
        "${trimgalore_dir}/${sample}_R1_001_val_1.fq.gz"
    do
        r2_candidate="${r1_candidate/_val_1.fq.gz/_val_2.fq.gz}"
        r2_candidate="${r2_candidate/_1_val_2.fq.gz/_2_val_2.fq.gz}"
        r2_candidate="${r2_candidate/_R1_val_2.fq.gz/_R2_val_2.fq.gz}"
        r2_candidate="${r2_candidate/_R1_001_val_2.fq.gz/_R2_001_val_2.fq.gz}"
        if [ -s "$r1_candidate" ] && [ -s "$r2_candidate" ]; then
            trimmed_done=true
            break
        fi
    done

    if [ "$trimmed_done" = true ]; then
        echo "跳过已完成修剪：$sample"
        continue
    fi

    echo "正在修剪样本：$sample"

    "$trim_galore" -q 25 --phred33 --length 50 \
        --paired --cores "$THREADS" --path_to_cutadapt "$cutadapt" \
        -o "$trimgalore_dir" "$fq1" "$fq2"
done

echo "✅ 修剪完成：共处理 $sample_count 个样本"

# -------------------------- 8. 修剪后FastQC --------------------------
echo -e "\n=== Step 4: 修剪后FastQC质控 ==="
trimmed_fastqs=( "$trimgalore_dir"/*val*.fq.gz )

if [ "${#trimmed_fastqs[@]}" -eq 0 ]; then
    echo "❌ 错误：未找到修剪后的 fq.gz 文件"
    exit 1
fi

"$fastqc" -t "$THREADS" -o "$fastqc_clean" "${trimmed_fastqs[@]}"
"$multiqc" -o "$fastqc_clean" "$fastqc_clean"
echo "✅ 修剪后质控完成"

# -------------------------- 9. BWA比对 --------------------------
echo -e "\n=== Step 5: BWA-MEM比对 ==="

tail -n +2 "$pair_tsv" | while IFS=$'\t' read -r sample fq1 fq2; do
    trimmed_fq1=""
    trimmed_fq2=""

    for candidate in \
        "${trimgalore_dir}/${sample}_1_val_1.fq.gz" \
        "${trimgalore_dir}/${sample}_R1_val_1.fq.gz" \
        "${trimgalore_dir}/${sample}_R1_001_val_1.fq.gz"
    do
        if [ -f "$candidate" ]; then
            trimmed_fq1="$candidate"
            break
        fi
    done

    for candidate in \
        "${trimgalore_dir}/${sample}_2_val_2.fq.gz" \
        "${trimgalore_dir}/${sample}_R2_val_2.fq.gz" \
        "${trimgalore_dir}/${sample}_R2_001_val_2.fq.gz"
    do
        if [ -f "$candidate" ]; then
            trimmed_fq2="$candidate"
            break
        fi
    done

    if [ -z "$trimmed_fq1" ] || [ -z "$trimmed_fq2" ]; then
        echo "❌ 错误：未找到修剪后的双端文件：$sample"
        echo "R1=$trimmed_fq1"
        echo "R2=$trimmed_fq2"
        exit 1
    fi

    sorted_bam="${sorted}${sample}.sorted.bam"
    log_out="${alignment}${sample}_bwa.log"

    if [ -s "$sorted_bam" ] && [ -s "${sorted_bam}.bai" ]; then
        echo "跳过已完成比对和排序：$sample"
        continue
    fi

    echo "正在比对并排序样本：$sample"

    "$bwa" mem -t "$THREADS" -M -T 30 \
        -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}_lib" \
        "$bwa_index" "$trimmed_fq1" "$trimmed_fq2" 2> "$log_out" | \
    "$samtools" sort -O bam -@ "$THREADS" -m "$SORT_MEM" -o "$sorted_bam" -

    "$samtools" index "$sorted_bam"
done

bwa_logs=( "$alignment"/*.log )
if [ "${#bwa_logs[@]}" -gt 0 ]; then
    "$multiqc" -o "$alignment" "$alignment"/*.log || true
fi
echo "✅ 比对完成"

# -------------------------- 10. 去重 --------------------------
echo -e "\n=== Step 7: Picard去重 ==="
cd "$sorted"

sorted_bams=( *.sorted.bam )
if [ "${#sorted_bams[@]}" -eq 0 ]; then
    echo "❌ 错误：未找到 sorted BAM 文件"
    exit 1
fi

for bam in "${sorted_bams[@]}"; do
    [ -f "$bam" ] || continue
    sample=$(basename "$bam" ".sorted.bam")
    dedup_bam="${dedup}${sample}.dedup.bam"
    dedup_metrics="${dedup}${sample}_dup_metrics.txt"

    if [ -s "$dedup_bam" ] && [ -s "${dedup_bam}.bai" ]; then
        echo "跳过已完成去重：$sample"
        continue
    fi

    "$picard" MarkDuplicates \
        I="$bam" \
        O="$dedup_bam" \
        M="$dedup_metrics" \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=SILENT \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

    "$samtools" index "$dedup_bam"
done

echo "✅ 去重完成"

# -------------------------- 11. 变异检测 --------------------------
echo -e "\n=== Step 8: 变异检测 ==="
cd "$dedup"

dedup_bams=( *.dedup.bam )
if [ "${#dedup_bams[@]}" -eq 0 ]; then
    echo "❌ 错误：未找到 dedup BAM 文件"
    exit 1
fi

for bam in "${dedup_bams[@]}"; do
    [ -f "$bam" ] || continue
    sample=$(basename "$bam" ".dedup.bam")
    raw_vcf_gz="${variant_raw}${sample}.raw.vcf.gz"
    mpileup_log="${variant_raw}${sample}_mpileup.log"

    if [ -s "$raw_vcf_gz" ] && { [ -s "${raw_vcf_gz}.csi" ] || [ -s "${raw_vcf_gz}.tbi" ]; }; then
        echo "跳过已完成变异检测：$sample"
        continue
    fi

    echo "正在检测变异：$sample"

    "$bcftools" mpileup \
        -a FORMAT/DP,FORMAT/AD \
        -f "$reference" \
        -d 1000 -q 20 -Q 20 -C 50 \
        "$bam" 2> "$mpileup_log" | \
    "$bcftools" call -mv -O z --ploidy 1 -o "$raw_vcf_gz"

    "$bcftools" index -f "$raw_vcf_gz"
done

echo "✅ 变异检测完成：$variant_raw"

# -------------------------- 13. 基础质量过滤 --------------------------
echo -e "\n=== Step 9: 基础质量过滤（保留多态突变） ==="
cd "$variant_raw"

raw_vcfs=( *.raw.vcf.gz )
if [ "${#raw_vcfs[@]}" -eq 0 ]; then
    echo "❌ 错误：未找到 raw VCF 文件"
    exit 1
fi

for vcf_gz in "${raw_vcfs[@]}"; do
    [ -f "$vcf_gz" ] || continue
    sample=$(basename "$vcf_gz" ".raw.vcf.gz")
    qc_vcf_gz="${variant_qc}${sample}.qc.vcf.gz"

    if [ -s "$qc_vcf_gz" ] && { [ -s "${qc_vcf_gz}.csi" ] || [ -s "${qc_vcf_gz}.tbi" ]; }; then
        echo "跳过已完成基础过滤：$sample"
        continue
    fi

    echo "正在基础过滤：$sample"

    "$bcftools" filter \
        -e "QUAL<${MIN_QUAL} || FORMAT/DP<${MIN_DP_CALL}" \
        -O z -o "$qc_vcf_gz" "$vcf_gz"

    if [ -s "$qc_vcf_gz" ]; then
        "$bcftools" index -f "$qc_vcf_gz"
    else
        echo "⚠️ 过滤后为空：$sample"
    fi
done

echo "✅ 基础质量过滤完成：$variant_qc"

# -------------------------- 14. 按AF拆分：polymorphic / highfreq --------------------------
echo -e "\n=== Step 10: 按AF拆分为多态突变和高频突变 ==="
cd "$variant_qc"

qc_vcfs=( *.qc.vcf.gz )
if [ "${#qc_vcfs[@]}" -eq 0 ]; then
    echo "❌ 错误：未找到 qc VCF 文件"
    exit 1
fi

for vcf_gz in "${qc_vcfs[@]}"; do
    [ -f "$vcf_gz" ] || continue
    sample=$(basename "$vcf_gz" ".qc.vcf.gz")

    poly_sites="${variant_tables}${sample}.polymorphic.sites.tsv"
    high_sites="${variant_tables}${sample}.highfreq.sites.tsv"
    poly_tsv="${variant_tables}${sample}.polymorphic.table.tsv"
    high_tsv="${variant_tables}${sample}.highfreq.table.tsv"

    poly_vcf_gz="${variant_poly}${sample}.polymorphic.vcf.gz"
    high_vcf_gz="${variant_high}${sample}.highfreq.vcf.gz"

    if [ -s "$poly_tsv" ] && [ -s "$high_tsv" ] \
        && { [ -s "$poly_vcf_gz" ] || [ ! -s "$poly_sites" ]; } \
        && { [ -s "$high_vcf_gz" ] || [ ! -s "$high_sites" ]; }; then
        echo "跳过已完成AF拆分：$sample"
        continue
    fi

    echo "正在AF拆分：$sample"

    "$python3_bin" "$split_py" \
        --vcf "$vcf_gz" \
        --poly_sites "$poly_sites" \
        --high_sites "$high_sites" \
        --poly_tsv "$poly_tsv" \
        --high_tsv "$high_tsv" \
        --min_dp "$MIN_DP_AF" \
        --min_af_poly "$MIN_AF_POLY" \
        --max_af_poly "$MAX_AF_POLY" \
        --min_af_high "$MIN_AF_HIGH"

    if [ -s "$poly_sites" ]; then
        "$bcftools" view -R "$poly_sites" -O z -o "$poly_vcf_gz" "$vcf_gz"
        "$bcftools" index -f "$poly_vcf_gz"
        echo "✅ polymorphic 结果已输出：$sample"
    else
        echo "⚠️ polymorphic 结果为空：$sample"
    fi

    if [ -s "$high_sites" ]; then
        "$bcftools" view -R "$high_sites" -O z -o "$high_vcf_gz" "$vcf_gz"
        "$bcftools" index -f "$high_vcf_gz"
        echo "✅ highfreq 结果已输出：$sample"
    else
        echo "⚠️ highfreq 结果为空：$sample"
    fi
done

echo "✅ AF拆分完成"

# -------------------------- 15. 注释函数 --------------------------
annotate_dir () {
    local input_dir="$1"
    local output_dir="$2"
    local label="$3"

    echo -e "\n=== Step 11 (${label}): snpEff注释 ==="
    cd "$input_dir"

    local vcfs=( *.vcf.gz )
    if [ "${#vcfs[@]}" -eq 0 ]; then
        echo "⚠️ ${label} 目录下没有 VCF 文件，跳过注释"
        return 0
    fi

    for vcf_gz in "${vcfs[@]}"; do
        [ -f "$vcf_gz" ] || continue
        sample=$(basename "$vcf_gz" ".vcf.gz")

        temp_vcf="${output_dir}${sample}.tmp.vcf"
        ann_vcf_gz="${output_dir}${sample}.annotated.vcf.gz"
        ann_log="${output_dir}${sample}_snpEff.log"
        stats_html="${output_dir}${sample}_snpEff_stats.html"

        if [ -s "$ann_vcf_gz" ] && { [ -s "${ann_vcf_gz}.csi" ] || [ -s "${ann_vcf_gz}.tbi" ]; }; then
            echo "跳过已完成注释：$sample"
            continue
        fi

        echo "正在注释：$sample"

        "$bcftools" view -O v "$vcf_gz" | \
        "$snpeff" ann -Xmx16g \
            "${snpeff_config_args[@]}" \
            -dataDir "$SNPEFF_DATA_DIR" \
            -stats "$stats_html" \
            -v "$CG_DB_NAME" - > "$temp_vcf" 2> "$ann_log"

        if [ -s "$temp_vcf" ]; then
            "$bcftools" view -O z -o "$ann_vcf_gz" "$temp_vcf"
            "$bcftools" index -f "$ann_vcf_gz"
            rm -f "$temp_vcf"
            echo "✅ 注释完成：$sample"
        else
            echo "⚠️ 注释结果为空：$sample"
            rm -f "$temp_vcf"
        fi
    done

    "$multiqc" -o "$output_dir" "$output_dir" || true
}

annotate_dir "$variant_poly" "$variant_poly_annot" "polymorphic"
annotate_dir "$variant_high" "$variant_high_annot" "highfreq"

# -------------------------- 16. 统计 --------------------------
echo -e "\n=== Step 12: 变异统计 ==="

for ann_dir in "$variant_poly_annot" "$variant_high_annot"; do
    cd "$ann_dir"

    ann_vcfs=( *.annotated.vcf.gz )
    if [ "${#ann_vcfs[@]}" -eq 0 ]; then
        echo "⚠️ $ann_dir 下没有 annotated VCF，跳过统计"
        continue
    fi

    for vcf_gz in "${ann_vcfs[@]}"; do
        [ -f "$vcf_gz" ] || continue
        sample=$(basename "$vcf_gz" ".annotated.vcf.gz")
        stats_out="${ann_dir}${sample}.variant_stats.txt"
        if [ -s "$stats_out" ]; then
            echo "跳过已完成统计：$sample"
            continue
        fi
        "$bcftools" stats "$vcf_gz" > "$stats_out"
    done

    "$multiqc" -o "$ann_dir" "$ann_dir" || true
done

# -------------------------- 17. 清理临时文件 --------------------------
echo -e "\n=== Step 13: 清理临时文件 ==="
rm -f "$pair_tsv"

# -------------------------- 18. 完成提示 --------------------------
echo -e "\n=============================================="
echo "✅ 谷氨酸棒杆菌重测序分析全流程完成（双输出版，多命名兼容，升级稳健版）！"
echo "=============================================="
echo "结果路径汇总："
echo "1. 原始质控：$fastqc_raw"
echo "2. 修剪后质控：$fastqc_clean"
echo "3. 清洁数据：$trimgalore_dir"
echo "4. 去重BAM：$dedup"
echo "5. 原始VCF：$variant_raw"
echo "6. 基础质量过滤VCF：$variant_qc"
echo "7. 多态突变VCF：$variant_poly"
echo "8. 高频突变VCF：$variant_high"
echo "9. 多态突变注释：$variant_poly_annot"
echo "10. 高频突变注释：$variant_high_annot"
echo "11. AF拆分表格：$variant_tables"
echo "12. 日志目录：$log_dir"
echo "=============================================="

exit 0
