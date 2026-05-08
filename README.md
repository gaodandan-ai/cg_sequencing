# Corynebacterium glutamicum Genome Resequencing Pipeline

谷氨酸棒杆菌（*Corynebacterium glutamicum*）全基因组重测序分析流程。项目包含从双端 FASTQ 质控、修剪、BWA 比对、BAM 排序去重、bcftools 变异检测，到按等位基因频率（AF）拆分 polymorphic / highfreq 突变并用 snpEff 注释的脚本。

## 目录结构

```text
.
├── 00_reference/              # 参考基因组 FASTA/GFF，可保留小型公共参考文件
├── 01_raw_data/               # 原始 FASTQ，本地放置，不上传 GitHub
├── 02_scripts/                # 主流程和辅助脚本
│   ├── CG_PrepareIndex.sh
│   ├── CG_GenomeReseq.sh
│   └── helpers/
├── 03_qc/                     # FastQC/MultiQC 输出，不上传
├── 04_clean_data/             # Trim Galore 输出，不上传
├── 05_alignment/              # SAM/BAM 输出，不上传
├── 06_variant_calling/        # VCF 和统计输出；默认仅保留 tables/*.tsv
├── 07_annotation/             # snpEff 自定义数据库，本地构建，不上传
└── 08_logs/                   # 运行日志，不上传
```

## 环境依赖

推荐使用 Conda/Mamba 创建环境：

```bash
mamba env create -f environment.yml
conda activate cg_seq
```

主要工具包括：BWA、Trim Galore、Cutadapt、FastQC、MultiQC、Samtools、Bcftools、snpEff、Picard 和 Python 3。

## 输入文件准备

1. 将参考基因组放入 `00_reference/`，默认文件名为：
   - `genome.fna`
   - `genome.gff`

2. 将双端测序数据放入 `01_raw_data/`。脚本支持以下命名：
   - `sample_1.fq.gz` / `sample_2.fq.gz`
   - `sample_1.fastq.gz` / `sample_2.fastq.gz`
   - `sample_R1.fq.gz` / `sample_R2.fq.gz`
   - `sample_R1_001.fq.gz` / `sample_R2_001.fq.gz`

3. 构建 snpEff 自定义数据库，并将数据库目录放在 `07_annotation/` 下。默认数据库名：

```text
Corynebacterium_glutamicum_ATCC13032
```

## 运行流程

先构建参考基因组索引：

```bash
bash 02_scripts/CG_PrepareIndex.sh
```

再运行重测序分析主流程：

```bash
bash 02_scripts/CG_GenomeReseq.sh
```

主流程会依次执行：

1. 原始 FASTQ 质控
2. Trim Galore 修剪
3. 修剪后质控
4. BWA-MEM 比对
5. SAM 转 BAM、排序和索引
6. Picard 去重
7. bcftools 变异检测
8. 基础质量过滤
9. 按 AF 拆分 polymorphic / highfreq 变异
10. snpEff 注释与 MultiQC 汇总

## 可配置参数

脚本默认使用项目内相对路径，也可以通过环境变量覆盖，例如：

```bash
REFERENCE=/path/to/genome.fna \
FASTQ_DIR=/path/to/fastq \
SNPEFF_DATA_DIR=/path/to/snpeff/data \
SNPEFF_CONFIG=/path/to/snpEff.config \
bash 02_scripts/CG_GenomeReseq.sh
```

工具路径也可覆盖，例如：

```bash
BWA=/path/to/bwa BCFTOOLS=/path/to/bcftools bash 02_scripts/CG_GenomeReseq.sh
```

关键阈值在 `02_scripts/CG_GenomeReseq.sh` 中设置：

```bash
MIN_DP_CALL=10
MIN_DP_AF=20
MIN_AF_POLY=0.01
MAX_AF_POLY=0.99
MIN_AF_HIGH=0.80
MIN_QUAL=30
```

## GitHub 上传说明

仓库已经通过 `.gitignore` 排除了原始 FASTQ、清洗 FASTQ、SAM/BAM、VCF、snpEff 二进制数据库、QC 报告和日志等大文件。默认可上传脚本、README、环境文件、参考 FASTA/GFF 以及 `06_variant_calling/tables/*.tsv` 这类小型结果表。

初始化并上传到 GitHub 的常用命令：

```bash
git init
git add .
git commit -m "Initial genome resequencing pipeline"
git branch -M main
git remote add origin git@github.com:<your-name>/<repo-name>.git
git push -u origin main
```

上传前建议检查将要提交的文件：

```bash
git status --short
git ls-files --others --ignored --exclude-standard
```
