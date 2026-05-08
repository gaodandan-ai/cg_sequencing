# 谷氨酸棒杆菌全基因组重测序分析流程

谷氨酸棒杆菌（*Corynebacterium glutamicum*）全基因组重测序分析流程。项目包含从双端 FASTQ 质控、修剪、BWA 比对、BAM 排序去重、bcftools 变异检测，到按等位基因频率（AF）拆分 polymorphic / highfreq 突变并用 snpEff 注释的脚本。

这个 README 尽量写成“边学边跑”的说明文档：先解释重测序分析的基本原理，再说明每一步脚本在做什么，最后给出常用终端命令和排错方法。

## 重测序分析原理

全基因组重测序（whole genome resequencing）通常不是从头组装一个新基因组，而是把样本测序得到的短序列 reads 比对到一个已知参考基因组上，然后寻找样本和参考基因组之间的差异。

可以把整个过程理解为：

```text
样本 DNA
  ↓ 测序
原始 FASTQ reads
  ↓ 质量控制和低质量序列修剪
clean FASTQ reads
  ↓ 比对到参考基因组
BAM 比对文件
  ↓ 变异检测
VCF 变异文件
  ↓ 过滤、分类、注释和统计
突变位点、突变频率、影响基因、氨基酸变化等结果
```

本项目主要回答这些问题：

- 每个样本相对于参考基因组有哪些 SNP 或 Indel？
- 每个位点的测序深度是多少？
- 支持突变的 reads 占比是多少，也就是 AF（allele frequency）？
- 哪些突变是高频突变，哪些突变属于多态突变？
- 突变位于哪些基因上，是否影响蛋白编码？

### 常见文件格式

`FASTQ`：测序仪输出的原始 reads 文件，常见后缀是 `.fq.gz` 或 `.fastq.gz`。文件中包含 read 序列和每个碱基的质量值。

`FASTA`：参考基因组序列文件。本项目默认使用 `00_reference/genome.fna`。

`GFF`：参考基因组注释文件，记录基因、CDS 等区域的位置。本项目默认使用 `00_reference/genome.gff`。

`BAM`：reads 比对到参考基因组后的压缩文件，记录每条 read 比到了哪里、比对质量如何等信息。

`VCF`：变异检测结果文件，记录突变位置、参考碱基、突变碱基、质量值、测序深度等信息。

### 本流程每一步在做什么

1. `FastQC` 检查原始 FASTQ 质量，例如碱基质量、接头污染、GC 含量和重复率。
2. `Trim Galore` / `Cutadapt` 去除接头和低质量 reads，得到更干净的 FASTQ。
3. 再次运行 `FastQC`，确认修剪后数据质量是否改善。
4. `BWA-MEM` 将 clean reads 比对到参考基因组。
5. `samtools sort` 直接把比对结果排序成 BAM，避免生成体积很大的 SAM 中间文件。
6. `Picard MarkDuplicates` 标记并去除 PCR 重复 reads。
7. `bcftools mpileup` 和 `bcftools call` 检测变异，生成 VCF。
8. `bcftools filter` 根据质量值和测序深度过滤低可信变异。
9. 辅助脚本 `split_vcf_by_af.py` 根据 AF 拆分 polymorphic 和 highfreq 变异。
10. `snpEff` 结合 GFF 注释文件解释突变影响，例如是否导致同义/非同义突变。

### AF 是什么

AF 是 allele frequency，即突变等位基因频率。可以简单理解为：

```text
AF = 支持突变的 reads 数 / 该位点总 reads 深度
```

例如某个位点总深度 `DP=100`，其中 `AO=85` 条 reads 支持突变，则：

```text
AF = 85 / 100 = 0.85
```

本流程默认：

- `polymorphic`：`0.01 <= AF < 0.99`
- `highfreq`：`AF >= 0.80`

注意：这两个分类默认不是互斥的。例如 `AF=0.85` 同时满足 polymorphic 和 highfreq。是否要把分类改成互斥，需要根据你的实验设计决定。

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

推荐使用 Conda激活环境：

```bash
conda env create -f environment.yml
conda activate cg_seq
```

激活环境后，可以用下面的命令检查常用工具是否能正常调用：

```bash
bwa
samtools --version
bcftools --version
fastqc --version
multiqc --version
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

检查输入文件是否放对位置：

```bash
ls -lh 00_reference/
ls -lh 01_raw_data/
ls -lh 07_annotation/
```

统计 FASTQ 文件数量：

```bash
find 01_raw_data -maxdepth 1 \( -name "*.fq.gz" -o -name "*.fastq.gz" \)
```

## 运行流程

进入项目目录：

```bash
cd ~/cg_sequencing
```

如果不确定自己在哪个目录：

```bash
pwd
ls
```

先构建参考基因组索引：

```bash
bash 02_scripts/CG_PrepareIndex.sh
```

构建完成后可以检查索引文件：

```bash
ls -lh 00_reference/index/
ls -lh 00_reference/genome.fna.fai
```

再运行重测序分析主流程：

```bash
bash 02_scripts/CG_GenomeReseq.sh
```

如果想把终端输出同时保存为日志，可以运行：

```bash
bash 02_scripts/CG_GenomeReseq.sh 2>&1 | tee 08_logs/sequencing_analysis/run.log
```

如果服务器 CPU 更多，可以临时指定线程数：

```bash
THREADS=8 SORT_MEM=2G bash 02_scripts/CG_GenomeReseq.sh
```

主流程会依次执行：

1. 原始 FASTQ 质控
2. Trim Galore 修剪
3. 修剪后质控
4. BWA-MEM 比对并直接生成排序 BAM
5. BAM 索引
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

也可以复制配置模板后按需修改：

```bash
cp config.env.example config.env
bash 02_scripts/CG_GenomeReseq.sh
```

关键阈值支持环境变量覆盖：

```bash
MIN_DP_CALL=10
MIN_DP_AF=20
MIN_AF_POLY=0.01
MAX_AF_POLY=0.99
MIN_AF_HIGH=0.80
MIN_QUAL=30
```

`THREADS` 和 `SORT_MEM` 可用于控制线程数和 `samtools sort` 单线程内存。

参数含义：

- `THREADS`：使用多少线程，机器资源充足时可以设置为 8、12 或 16。
- `SORT_MEM`：`samtools sort` 每个线程可用内存，例如 `1G` 或 `2G`。
- `REFERENCE`：参考基因组 FASTA 文件路径。
- `BWA_INDEX`：BWA 索引前缀，默认是 `00_reference/index/genome`。
- `FASTQ_DIR`：原始 FASTQ 文件目录。
- `SNPEFF_DATA_DIR`：snpEff 数据库目录。
- `CG_DB_NAME`：snpEff 数据库名称。
- `MIN_DP_CALL`：变异检测基础过滤使用的最低测序深度。
- `MIN_DP_AF`：计算和拆分 AF 时使用的最低测序深度。
- `MIN_AF_POLY` / `MAX_AF_POLY`：polymorphic 变异的 AF 范围。
- `MIN_AF_HIGH`：highfreq 变异的最低 AF。
- `MIN_QUAL`：VCF 变异质量值过滤阈值。

## 脚本特性

- 主流程直接调用 `02_scripts/helpers/` 中的辅助 Python 脚本，不会在运行时覆盖源码。
- BWA 比对结果直接通过管道进入 `samtools sort`，不再生成体积很大的 SAM 中间文件。
- 主要步骤支持断点续跑：已存在且完成索引的 BAM/VCF/注释/统计结果会自动跳过。

## 如何查看结果

### 质控报告

原始数据质控结果：

```bash
ls 03_qc/raw/
```

修剪后质控结果：

```bash
ls 03_qc/trimmed/
```

重点查看 `multiqc_report.html`。如果在服务器终端里看不到网页，可以把 HTML 下载到本地电脑后用浏览器打开。

### 清洗后的 FASTQ

```bash
ls -lh 04_clean_data/
```

### BAM 比对结果

排序后的 BAM：

```bash
ls -lh 05_alignment/sorted_bam/
```

去重后的 BAM：

```bash
ls -lh 05_alignment/dedup_bam/
```

查看某个样本 BAM 的比对统计：

```bash
samtools flagstat 05_alignment/dedup_bam/sample.dedup.bam
```

把 `sample` 替换成真实样本名，例如：

```bash
samtools flagstat 05_alignment/dedup_bam/C0.dedup.bam
```

### VCF 变异结果

原始 VCF：

```bash
ls -lh 06_variant_calling/raw/
```

过滤后 VCF：

```bash
ls -lh 06_variant_calling/quality_filtered/
```

polymorphic 和 highfreq VCF：

```bash
ls -lh 06_variant_calling/polymorphic/
ls -lh 06_variant_calling/highfreq/
```

### AF 拆分表格

```bash
ls -lh 06_variant_calling/tables/
```

查看前几行：

```bash
head 06_variant_calling/tables/C0.highfreq.table.tsv
head 06_variant_calling/tables/C0.polymorphic.table.tsv
```

表格字段含义：

```text
CHROM  染色体或 contig
POS    位置
REF    参考碱基
ALT    突变碱基
DP     总测序深度
AO     支持突变的 reads 数
AF     突变等位基因频率
```

### snpEff 注释结果

```bash
ls -lh 06_variant_calling/annotated/highfreq/
ls -lh 06_variant_calling/annotated/polymorphic/
```

常见输出：

```text
sample.highfreq.annotated.vcf.gz
sample.polymorphic.annotated.vcf.gz
sample_snpEff_stats.html
sample.variant_stats.txt
```


## 常见问题

### 找不到 FASTQ 文件

报错可能类似：

```text
错误：01_raw_data 下没有 .fastq.gz 或 .fq.gz 文件
```

检查：

```bash
ls -lh 01_raw_data/
```

确认文件后缀是 `.fq.gz` 或 `.fastq.gz`，不要解压成 `.fq`。

### 双端文件不完整

脚本要求 R1 和 R2 成对出现。推荐命名：

```text
sample_1.fq.gz
sample_2.fq.gz
```

不要混用：

```text
sample_1.fq.gz
sample_R2.fq.gz
```

### 找不到软件

先激活环境：

```bash
conda activate cg_seq
```

再检查：

```bash
which bwa
which samtools
which bcftools
```

### snpEff 数据库不存在

检查数据库目录：

```bash
ls -lh 07_annotation/
```

确认 `config.env` 中的 `CG_DB_NAME` 和实际目录名一致。

### 中途失败后是否需要从头运行

不需要。脚本支持断点续跑，修复错误后重新运行主脚本即可：

```bash
bash 02_scripts/CG_GenomeReseq.sh
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

