# nextflow-metagenomics

[English](README.md) | **中文**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Nextflow 26.04.4](https://img.shields.io/badge/Nextflow-26.04.4-0dc09d?logo=nextflow&logoColor=white)](https://www.nextflow.io/)
![conda](https://img.shields.io/badge/conda-supported-3EB049?logo=anaconda&logoColor=white)
![Docker](https://img.shields.io/badge/Docker-supported-2496ED?logo=docker&logoColor=white)
![Singularity](https://img.shields.io/badge/Singularity-supported-5E97F6)

基于 [Nextflow DSL2](https://www.nextflow.io/) 的模块化、可复现鸟枪法宏基因组学
工作流：物种分类、功能谱分析、跨样本整合、MAG 重建与结果可视化。

**项目状态：** 完整流程已实现（输入 → 预处理 → read-based ∥ 组装 → MAG →
质控 → 去冗余 → 分类 → 基因预测 → 注释 → 丰度 → 整合 → 报告 → 可视化）。

> **验证范围**（哪些已真正执行、哪些待验证）：预处理、组装（MEGAHIT）、
> 比对、分箱、基因预测、丰度、整合、MultiQC 已在合成测试数据上真实运行验证。
> Kraken2、Bracken、HUMAnN、CheckM2、dRep、GTDB-Tk、DIAMOND、eggNOG-mapper、
> RGI 依赖数 GB 级数据库，开发机未提供 —— 仅通过 **stub-run** 验证（通道拓扑、
> 参数守卫、输出结构）。跨样本合并逻辑与绘图（除无库依赖的
> CoverM 丰度热图外）同样 database-dependent，由 stub/合成表验证。**真实数据库
> 执行与容器实机运行均未验证。** 详见 [Limitations](#limitations)。

---

## 目录

- [项目简介](#项目简介)
- [Workflow](#workflow)
- [仓库结构](#仓库结构)
- [核心特性](#核心特性)
- [快速开始](#快速开始)
- [输入与输出](#输入与输出)
- [配置与数据库](#配置与数据库)
- [测试与可重复性](#测试与可重复性)
- [Limitations](#limitations)
- [工具与引用](#工具与引用)
- [文档](#文档)
- [License](#license)

## 项目简介

`nextflow-metagenomics` 是一个用于 shotgun metagenomics 分析的 Nextflow DSL2
工作流。它消费原始 FASTQ，执行 QC 与可选宿主去除，随后从同一份 clean reads
**并行**运行 read-based 分析（物种分类、丰度估计、功能谱分析）与组装/MAG
路径（组装、MAG 重建、质控、去冗余、分类、基因预测、功能注释、丰度）。逐样本
read-based 结果合并为跨样本矩阵，全部结果 join 为核心表，最终生成 MultiQC
报告与一组可视化图。

项目按分层、单一职责的 Nextflow DSL2 结构组织：薄入口 → 一个编排 workflow →
每阶段一个子工作流 → 本地 process 模块，`tuple(meta, …)` 元数据贯穿从输入
到报告的每条通道。

## Workflow

| 阶段 | 工具 |
| ---- | ---- |
| 质控与宿主去除 | FastQC → fastp → Bowtie2 + samtools |
| Read-based 分析 | Kraken2 → Bracken ∥ HUMAnN |
| 组装 | MEGAHIT ∥ metaSPAdes → QUAST |
| 比对与覆盖度 | Bowtie2 → samtools → depth |
| 分箱 | MetaBAT2 |
| MAG 质控 | CheckM2 |
| 去冗余 | dRep |
| 物种分类 | GTDB-Tk |
| 基因预测 | Prodigal |
| 功能注释 | DIAMOND ∥ eggNOG-mapper ∥ RGI (CARD) |
| MAG 丰度 | CoverM |
| 整合 | 自研 Python —— 两张核心表 |
| 报告 | MultiQC + 7 张图（自研 Python） |

![nextflow-metagenomics workflow overview](docs/images/workflow_compact.png)

- **并行：** read-based 分析与组装/MAG 链从 clean reads 分叉，彼此无依赖。
- **串联：** 组装 → 比对 → 分箱，以及 MAG 质控 → 去冗余 → 分类（dRep 打分消费
  QC 输出，GTDB-Tk 只分类去冗余后的代表集）。
- **扇出：** 基因预测 / 注释 / 丰度并行消费代表集；后续整合阶段汇聚其汇总表。
- 完整竖版见 [docs/images/workflow.png](docs/images/workflow.png)。

## 仓库结构

```text
nextflow-metagenomics/
├── main.nf                    # 薄入口 → workflows/mag.nf
├── workflows/mag.nf           # 主 workflow：全部阶段在此接线
├── nextflow.config            # 参数定义 + profiles + db_dir 派生
├── conf/                      # base / conda / docker / singularity / slurm / test
├── subworkflows/local/        # 14 个子工作流（每阶段一个）
├── modules/local/             # 45 个本地模块 / 47 个 process（单一职责）
├── bin/                       # 14 个 Python 解析/编排脚本 + 1 个绘图工具脚本
├── assets/                    # 哨兵文件（如 empty.tsv）
├── test/                      # 合成数据、测试配置、单元测试
├── docs/                      # architecture / channels / workflow / input / output / database / parameters
├── environment.yml            # 钉版开发环境（nf-meta）
├── setup_env.sh               # 建环境 + 安装后修复
└── LICENSE
```

## 核心特性

- **Nextflow DSL2 分层架构** —— 薄 `main.nf` → `workflows/mag.nf` → 14 个子工作流
  → 45 个本地模块（47 个 process 定义）。每个 process 单一职责（`FASTQC` /
  `FASTP`，而非 `RUN_ALL_ANALYSIS`），并各自声明 conda/container 指令、
  publishDir、stub 块与 `versions.yml`。

- **Channel 与 metadata 契约** —— 每条通道都是 `tuple(meta, …)`；
  `meta = [id, group, batch, host, single_end]` 贯穿并随阶段扩展（组装后追加
  `assembler`、`assembly_mode`、`samples`）。

- **显式数据流设计（`join` / `combine` / `collect` / `groupTuple`）：**
  - `combine` 把 Bracken 按分类层级 scatter（样本 × 层级）。
  - `combine(by: 0)` 以样本身份为键把 contigs 与「组装出它的那些 reads」配对
    （在 `--assembler both` 下依然正确）。
  - `groupTuple` + `groupKey` 按组装单元聚合排序 BAM。
  - `collectFile` / `toSortedList` 聚合批级输入并保证确定性排序，使 `-resume`
    命中缓存。
  - `flatMap` 把一个 bin 目录展开为每个 MAG 一条记录；空通道是合法状态
    （0-bin 样本）。

- **参数化与 skip 策略：**
  - 约 80 个参数全部定义于 `nextflow.config`；数据库/索引路径只经 `params.*`
    传入（无硬编码路径）。
  - 每个工具有 `--skip_*` 开关、每个分支有总开关；read-based 数据库缺失时
    **告警并跳过**该并行分支，MAG 级数据库（CheckM2 / GTDB-Tk / 注释）缺失时
    **显式报错** —— 不静默降级、不虚构产物。

- **资源管理** —— `conf/base.config` 定义四档 label
  （`process_single/low/medium/high`），由 `resourceLimits` 封顶；对高内存工具
  （KRAKEN2、GTDBTK 各 64 GB）作逐工具覆盖，并对资源类退出码做带退避的重试。

- **可重复性与运行 profile** —— 版本钉死（`environment.yml` + 每个 process 的
  conda/container 指令）、固定 MetaBAT2 随机种子、每个 process 产出
  `versions.yml`、支持 `-stub-run`。提供 `conda` / `docker` / `singularity` /
  `slurm` / `test` 等 profile，以及以 `--db_dir` 一棵树派生各数据库路径的集中
  布局。

- **逐样本 → 跨样本整合** —— 把逐样本 Bracken/HUMAnN 表合并为样本×taxa 与
  样本×pathway 矩阵（另含 Bray–Curtis 距离矩阵）；再把 QC/分类/注释/丰度 join
  成两张核心表，并作为只读消费层渲染 7 张图。

## 快速开始

安装 —— 请用 `setup_env.sh`（它含 `mamba env create` 单独无法表达的安装后修复；
细节见 `environment.yml` 顶部注释）：

```bash
bash setup_env.sh          # 创建 + 修复 + 逐工具启动验证
conda activate nf-meta
```

运行：

```bash
# 标准运行 (read-based 与组装分支; 未提供数据库的分支自动跳过并告警)
nextflow run main.nf --input samplesheet.csv

# 指定 batch id (输出到 results/<batch_id>/)
nextflow run main.nf --input samplesheet.csv --batch_id batch_001

# 测试运行: 内置合成数据 (需先建好 nf-meta 环境)
nextflow run main.nf -profile test \
  --skip_read_based --skip_mag_qc --skip_dereplication \
  --skip_taxonomy --skip_annotation

# Stub 运行: 只验证通道拓扑与输出结构, 无需数据或数据库
# (dummy/fake DB 参数让每条分支都被调度, 但不跑真实工具)
nextflow run main.nf -profile test -stub-run \
  --kraken2_db test/data/dummy_dbs/kraken2 --bracken_db test/data/dummy_dbs/kraken2 \
  --humann_db test/data/dummy_dbs/humann \
  --metaphlan_db test/data/dummy_dbs/humann/metaphlan \
  --checkm2_db /fake --gtdbtk_db /fake --diamond_db /fake \
  --eggnog_db /fake --card_db /fake

# 其它运行时 / 调度器 (按需补充数据库)
nextflow run main.nf --input samplesheet.csv -profile docker
nextflow run main.nf --input samplesheet.csv -profile singularity     # Apptainer/Singularity
nextflow run main.nf --input samplesheet.csv -profile slurm
nextflow run main.nf --input samplesheet.csv -profile conda

# 数据库集中布局: 一条参数代替八条
nextflow run main.nf --input samplesheet.csv --db_dir /path/to/databases
```

Profile 可叠加（如 `-profile test,docker -stub-run`）。

## 输入与输出

### 输入

`--input` 指向 samplesheet CSV（每行一个样本；`fastq_2` 留空即单端）：

```csv
sample,fastq_1,fastq_2,group,batch,host
S01,/data/reads/S01_R1.fastq.gz,/data/reads/S01_R2.fastq.gz,case,batch01,human
S02,/data/reads/S02_R1.fastq.gz,/data/reads/S02_R2.fastq.gz,control,batch01,human
```

| 列 | 含义 |
| ---- | ---- |
| `sample` | 唯一样本 id（即 `meta.id`） |
| `fastq_1` | 正向 reads（R1）路径，必填 |
| `fastq_2` | 反向 reads（R2）路径；留空即单端 |
| `group` | 样本分组标签（如 `case` / `control`） |
| `batch` | 批次标签，随 metadata 传递 |
| `host` | 宿主基因组标签，用于宿主去除（如 `human`） |

`single_end` 会自动推导（`fastq_2` 留空即单端）并写入验证后的输出。
`CHECK_SAMPLESHEET` process 验证输入、把相对路径补全为绝对路径，并产出
`00_metadata/validated_samplesheet.csv`。详见 [docs/input.md](docs/input.md)。

### 输出

结果按 batch 组织在 `results/<batch_id>/` 下的固定编号目录（`00_metadata/` …
`14_integrated/`、`99_multiqc/`）。两张核心交付：

- `mag_metadata.tsv` —— QC / 分类 / 基因组统计 / 丰度，行 = 代表 MAG。
- `mag_functional_annotation.tsv` —— KO / COG / GO / Pathway / ARG，行 = gene。

跨样本整合在 `03_taxonomy/combined/` 与 `04_function/combined/` 下新增跨样本
矩阵；绘图在相关阶段目录下新增 `figures/` 子目录（7 张 PNG）。逐文件
清单见 [docs/output.md](docs/output.md)。

## 配置与数据库

### 配置与 Skip 策略

全部参数定义于 `nextflow.config`；完整清单（按阶段分组、含默认值）见
[docs/parameters.md](docs/parameters.md)。要点：

- **publishDir 闭包**延迟求值 `${params.outdir}/${params.batch_id}/...`，使
  profile/CLI 覆盖生效。
- **skip 语义具传递性**：跳过上游会使下游以空通道收尾（不报错）。
  `--skip_mag_qc` 必须连带 `--skip_dereplication`（dRep 打分消费 QC 输出）。
- **数据库**经 `--*_db` 传入（或由单个 `--db_dir` 树派生）；见
  [docs/database.md](docs/database.md)。

### 数据库

| 分支 / 阶段 | 数据库参数 | 缺失时行为 |
| --- | --- | --- |
| Kraken2 / Bracken | `--kraken2_db`（`--bracken_db` 默认复用它） | 告警并跳过该分支 |
| HUMAnN | `--humann_db`（+ `--metaphlan_db`） | 告警并跳过该分支 |
| CheckM2 | `--checkm2_db`（~3 GB） | 报错（或 `--skip_mag_qc`） |
| GTDB-Tk | `--gtdbtk_db`（R220+，~110 GB） | 报错（或 `--skip_taxonomy`） |
| DIAMOND | `--diamond_db`（NR `.dmnd`） | 报错（或 `--skip_diamond`） |
| eggNOG-mapper | `--eggnog_db`（~40+ GB） | 报错（或 `--skip_eggnog`） |
| RGI / CARD | `--card_db`（card.json） | 报错（或 `--skip_rgi`） |
| Pathway（可选） | `--pathway_db`（KO→pathway） | 列留空 |

dRep、Prodigal、CoverM 无外部数据库依赖。完整清单（体量、获取方式、
`--db_dir` 布局）见 [docs/database.md](docs/database.md)。

## 测试与可重复性

### 测试与验证

两条主要回归命令（逐阶段矩阵见 [docs/workflow.md](docs/workflow.md) 与
`STATUS.md`）：

```bash
# Stub 全流程 — 62/62 tasks, 无真实工具执行, 无数据库依赖;
# 验证通道拓扑、参数守卫与输出结构
nextflow run main.nf -profile test -stub-run --batch_id demo_stub \
  --kraken2_db test/data/dummy_dbs/kraken2 --bracken_db test/data/dummy_dbs/kraken2 \
  --humann_db test/data/dummy_dbs/humann \
  --metaphlan_db test/data/dummy_dbs/humann/metaphlan \
  --checkm2_db /fake --gtdbtk_db /fake --diamond_db /fake \
  --eggnog_db /fake --card_db /fake

# 真实全流程 — 34/34 tasks, 需 nf-meta 环境; read-based 与数据库依赖阶段
# 被 skip (数据库本地未提供)
nextflow run main.nf -profile test --batch_id demo_real \
  --skip_read_based --skip_mag_qc --skip_dereplication \
  --skip_taxonomy --skip_annotation
```

- 测试数据为 ~46 kb 合成宏基因组（2 个 PE 样本，经
  `test/data/make_test_data.py` 可复现）。
- 解析/编排脚本有单元测试覆盖（`test/test_merge_read_based.py`、
  `test/test_plot_results.py`，以及 `STATUS.md` 中记录的各脚本套件）。
- `-stub-run` 运行每个 process 的 stub（解析类 process 仍跑真实解析脚本），
  因此无需工具或数据库即可检验通道接线与输出 schema。

### 可重复性与部署

- **版本钉死：** `environment.yml` 钉住每个工具（Nextflow 26.04.4 —— 全部验证
  实际运行的版本）；每个 process 还声明 conda 钉版与 container tag。
- **`-resume` 幂等：** 批级聚合用 `collectFile(sort)` + `toSortedList` 保证
  task hash 确定（stub 与真实运行均实测缓存命中）。
- **固定种子：** `metabat2_seed=42` 使 MAG ID 稳定。
- **`versions.yml`** 每个 process 产出，MultiQC 汇总。
- **无硬编码路径：** 所有数据库/索引路径来自参数；gitignored 的
  `conf/local.config` 是唯一存放个人本地路径的位置。
- **Profile：** 本地（无 profile，工具取自 PATH）/ `conda` / `docker` /
  `singularity` / `slurm` / `test` —— 见 [conf/](conf/) 与
  [docs/architecture.md](docs/architecture.md)。

## Limitations

- **database-dependent 分支真实运行未验证：** read-based（Kraken2 标准库数十
  GB）、CheckM2（~3 GB）、dRep、GTDB-Tk（~110 GB）与注释三支（NR / eggNOG /
  CARD）仅 stub-run 验证，未虚构任何结果。
- **跨样本合并逻辑与绘图**（除 CoverM 丰度热图外）由 stub/合成表
  验证，待真实数据库。
- **metaSPAdes 真实运行未验证**（MEGAHIT 已真实验证）。
- **容器路径未实机验证：** 开发机无 Docker/Apptainer/Singularity/sbatch；已核验
  镜像 tag 存在性与 profile 叠加语法，实际拉取/容器内运行待有引擎环境。
- **测试规模：** 46 kb 合成数据、2 样本、1 个 bin —— 验证的是流程正确性与数值
  交叉一致，而非真实宏基因组的生物学结论。
- **尚未实现：** `assembly_mode = coassembly`（显式报错；通道已预留
  `meta.assembly_mode` / `meta.samples`）以及 MetaBAT2 之外的 binner
  （`binner` 参数为 V2 预留）。

## 工具与引用

本工作流编排了多个第三方工具；在你自己的工作中使用结果时，请引用其原始文献。
本工作流未提出新方法。

| 工具 | 作用 | 主要文献 |
| ---- | ---- | -------- |
| [FastQC](https://github.com/s-andrews/FastQC) | read 质控 | Andrews, 2010 |
| [fastp](https://github.com/OpenGene/fastp) | read 去接头 / 质控 | Chen 等, 2018, Bioinformatics |
| [Bowtie2](https://github.com/BenLangmead/bowtie2) | 宿主去除 + 比对 | Langmead & Salzberg, 2012, Nat. Methods |
| [samtools](https://github.com/samtools/samtools) | BAM 处理 | Danecek 等, 2021, GigaScience |
| [Kraken2](https://github.com/DerrickWood/kraken2) | 物种分类 | Wood, Lu & Langmead, 2019, Genome Biol. |
| [Bracken](https://github.com/jenniferlu717/Bracken) | 丰度估计 | Lu 等, 2017, PeerJ CS |
| [HUMAnN](https://github.com/biobakery/humann) | 功能谱分析 | Franzosa 等, 2018, Nat. Methods |
| [MEGAHIT](https://github.com/voutcn/megahit) | 组装 | Li 等, 2015, Bioinformatics |
| [metaSPAdes](https://github.com/ablab/spades) | 组装 | Nurk 等, 2017, Genome Res. |
| [QUAST](https://github.com/ablab/quast) | 组装质控 | Gurevich 等, 2013, Bioinformatics |
| [MetaBAT2](https://bitbucket.org/berkeleylab/metabat) | 分箱 + 深度 | Kang 等, 2019, PeerJ |
| [CheckM2](https://github.com/chklovski/CheckM2) | MAG 质控 | Chklovski 等, 2023, Nat. Methods |
| [dRep](https://github.com/MrOlm/drep) | 去冗余 | Olm 等, 2017, ISME J |
| [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) | MAG 物种分类 | Chaumeil 等, 2022, Bioinformatics |
| [Prodigal](https://github.com/hyattpd/Prodigal) | 基因预测 | Hyatt 等, 2010, BMC Bioinformatics |
| [DIAMOND](https://github.com/bbuchfink/diamond) | 蛋白比对 | Buchfink 等, 2021, Nat. Methods |
| [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) | 功能注释 | Cantalapiedra 等, 2021, Mol. Biol. Evol. |
| [RGI](https://github.com/arpcard/rgi) / CARD | 耐药基因 | Alcock 等, 2023, Nucleic Acids Res. |
| [CoverM](https://github.com/wwood/CoverM) | MAG 丰度 | Aroney 等, 2024, Bioinformatics |
| [MultiQC](https://github.com/MultiQC/MultiQC) | 报告汇总 | Ewels 等, 2016, Bioinformatics |

## 文档

| 文档 | 内容 |
| --- | --- |
| [docs/architecture.md](docs/architecture.md) | 分层架构、设计决策、工作流 DAG |
| [docs/workflow.md](docs/workflow.md) | 各阶段流程、mermaid 图、skip 交互语义 |
| [docs/channels.md](docs/channels.md) | 通道契约与操作符模式 |
| [docs/input.md](docs/input.md) | 输入格式与验证 |
| [docs/output.md](docs/output.md) | 输出目录逐文件清单 |
| [docs/database.md](docs/database.md) | 数据库清单、体量、`--db_dir` 布局 |
| [docs/parameters.md](docs/parameters.md) | 完整参数参考（按阶段分组） |

## License

本项目基于 MIT License 发布，仅用于学习与教学用途，详见 [LICENSE](LICENSE)。
