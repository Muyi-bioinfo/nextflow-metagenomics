# 参数参考

本文档是完整参数清单（按 Phase 分组）。**权威来源是 `nextflow.config` 的
`params` 块**（每行均带注释）—— 本文档只是便于阅读的整理，二者冲突时以
`nextflow.config` 为准。`—` 表示 `null` / 空。

## 通用

| 参数 | 默认值 | 用途 |
|------|--------|------|
| `input` | `—` | samplesheet CSV 路径（必填） |
| `outdir` | `results` | 结果根目录 |
| `batch_id` | `batch_<date>_001` | batch 标识，输出到 `results/<batch_id>/` |
| `threads` | `4` | 默认线程数 |
| `max_cpus` / `max_memory` / `max_time` | `16` / `10.GB` / `24.h` | 资源硬上限（resourceLimits） |

## 预处理（Phase 3）

| 参数 | 默认值 | 用途 |
|------|--------|------|
| `host_index` | `—` | Bowtie2 宿主索引前缀；缺失时跳过宿主去除 |
| `skip_host_removal` | `false` | 显式跳过宿主去除 |
| `fastp_qualified_quality` | `15` | 合格碱基质量阈值 |
| `fastp_unqualified_percent` | `40` | 允许的不合格碱基百分比上限 |
| `fastp_min_length` | `50` | 最短 read 长度 |
| `fastp_cut_mean_quality` | `20` | 滑窗平均质量修剪阈值 |
| `fastp_dedup` | `false` | fastp 去重 |
| `save_trimmed` | `true` | 发布 fastp 清洁 reads |
| `save_host_removed` | `false` | 发布去宿主 reads（体积大） |
| `skip_fastqc` | `false` | 跳过 FastQC |
| `skip_multiqc` | `false` | 跳过 MultiQC 报告 |

## Read-based（Phase 4）

| 参数 | 默认值 | 用途 |
|------|--------|------|
| `skip_read_based` | `false` | 总开关：跳过整个 read-based 分支 |
| `skip_kraken2` / `skip_bracken` / `skip_humann` | `false` | 逐工具开关 |
| `kraken2_confidence` | `0.0` | 置信度阈值 0–1 |
| `kraken2_min_base_quality` | `0` | 参与 k-mer 匹配的最低碱基质量 |
| `kraken2_min_hit_groups` | `2` | 判定已分类所需的最少命中组数 |
| `kraken2_memory_mapping` | `false` | true 时库不载入内存（省内存、更慢） |
| `save_kraken2_output` | `false` | 发布逐 read 分类结果（可达数 GB） |
| `bracken_levels` | `S,G` | 丰度估计层级（逗号分隔） |
| `bracken_read_length` | `100` | 必须与 bracken-build 时的读长一致 |
| `bracken_threshold` | `10` | 低于此 read 数的分类单元不重分配 |
| `humann_nucleotide_db` / `humann_protein_db` / `metaphlan_db` | `—` | HUMAnN 三个库（可由 `humann_db` 父目录推导；MetaPhlAn 必须显式给出） |
| `humann_args` | `''` | 追加给 humann 的参数 |

## 组装（Phase 5）

| 参数 | 默认值 | 用途 |
|------|--------|------|
| `skip_assembly` | `false` | 总开关：跳过组装分支 |
| `skip_quast` | `false` | 跳过组装 QC（连带跳过 assembly_summary.tsv） |
| `assembler` | `megahit` | `megahit` \| `metaspades` \| `both`（可选替代组装器，非两级流水线） |
| `assembly_mode` | `single` | `coassembly` 尚未实现（显式报错） |
| `megahit_min_contig_len` | `200` | 输出 contig 最短长度 |
| `megahit_k_list` / `megahit_preset` | `—` | k-mer 列表 / 预设（互斥） |
| `megahit_min_count` | `—` | 最小 (k+1)-mer 丰度（默认 2；低深度可设 1） |
| `megahit_args` | `''` | 追加参数 |
| `metaspades_k` | `—` | k-mer 列表（如 `21,33,55`） |
| `metaspades_args` | `''` | 追加参数 |
| `save_assembly_graph` | `false` | 发布 metaSPAdes GFA（大） |
| `quast_min_contig` | `500` | 计入统计的最短 contig 长度 |
| `quast_args` | `''` | 追加参数 |

## 比对与覆盖度（Phase 6）

| 参数 | 默认值 | 用途 |
|------|--------|------|
| `skip_mapping` | `false` | 跳过比对（连带分箱无输入） |
| `bowtie2_build_args` / `bowtie2_map_args` | `''` | 追加参数 |
| `save_bowtie2_index` | `false` | 发布 contigs 的 Bowtie2 索引 |
| `save_bam` | `false` | 发布排序 BAM+BAI（可达数十 GB） |
| `coverage_min_contig_len` / `coverage_min_depth` | `—` | jgi 覆盖度过滤（默认沿用工具默认值；长度过滤是 MetaBAT2 的职责） |

## 分箱（Phase 7）

| 参数 | 默认值 | 用途 |
|------|--------|------|
| `skip_binning` | `false` | 跳过分箱（连带 MAG 分析无输入） |
| `binner` | `metabat2` | V2 预留 maxbin2/concoct/dastools |
| `metabat2_min_contig_len` | `—` | `-m`（MetaBAT2 默认 2500） |
| `metabat2_min_bin_size` | `—` | `-s`（默认 200 kb） |
| `metabat2_seed` | `42` | 固定随机种子（MAG ID 可重现） |
| `metabat2_args` | `''` | 追加参数 |
| `save_unbinned` | `false` | 发布 .unbinned.fa |

## MAG QC / 去冗余 / 分类（Phase 8–10）

| 参数 | 默认值 | 用途 |
|------|--------|------|
| `skip_mag_qc` | `false` | 跳过 CheckM2（**必须连带 `--skip_dereplication`**） |
| `checkm2_args` | `''` | 追加参数 |
| `mag_min_completeness` | `50` | 完整度阈值（%），过滤 qualified MAG |
| `mag_max_contamination` | `10` | 污染度阈值（%） |
| `skip_dereplication` | `false` | 跳过 dRep：代表集 = 全部 qualified MAG，成员表为恒等映射 |
| `drep_args` | `''` | 追加参数（勿覆盖 -g/--genomeInfo/-p） |
| `skip_taxonomy` | `false` | 跳过 GTDB-Tk |
| `gtdbtk_args` | `''` | 追加参数（勿覆盖 --genome_dir/--out_dir/-x/--prefix） |

## 基因预测 / 注释（Phase 11–12）

| 参数 | 默认值 | 用途 |
|------|--------|------|
| `skip_gene_prediction` | `false` | 跳过 Prodigal（连带注释无输入） |
| `prodigal_args` | `''` | 追加参数（勿覆盖 -i/-a/-d/-f/-o/-p） |
| `skip_annotation` | `false` | 总开关：跳过功能注释 |
| `skip_diamond` / `skip_eggnog` / `skip_rgi` | `false` | 逐工具开关 |
| `diamond_args` / `eggnog_args` / `rgi_args` | `''` | 追加参数 |

## 丰度 / 整合（Phase 13–14）

| 参数 | 默认值 | 用途 |
|------|--------|------|
| `skip_abundance` | `false` | 跳过 CoverM |
| `coverm_method` | `relative_abundance` | CoverM `--methods`（保持单方法；Phase 14 不支持多方法矩阵） |
| `coverm_args` | `''` | 追加参数 |
| `skip_integration` | `false` | 跳过结果整合 |
| `pathway_db` | `—` | KO→pathway 两列映射（可选；缺失时 Pathway 列留空） |

## 可视化（Phase 21）

| 参数 | 默认值 | 用途 |
|------|--------|------|
| `plot_pathway_top` | `50` | 通路丰度热图的 top-N 通路数 |

## 数据库参数

| 参数 | 用途 | 缺失时行为 |
|------|------|-----------|
| `db_dir` | 集中布局根目录（见 docs/database.md） | `—` |
| `kraken2_db` / `bracken_db` | Kraken2 分类 / Bracken 丰度（默认复用 kraken2 库） | **告警并跳过**该分支 |
| `humann_db` | HUMAnN 功能谱父目录（含 chocophlan/uniref；metaphlan 单独传） | **告警并跳过**该分支 |
| `checkm2_db` | CheckM2 MAG QC（~3 GB） | **报错**（或 `--skip_mag_qc`） |
| `gtdbtk_db` | GTDB-Tk 分类（R220+，解压约 110 GB） | **报错**（或 `--skip_taxonomy`） |
| `diamond_db` | DIAMOND blastp（NR `.dmnd`，数十 GB） | **报错**（或 `--skip_diamond`） |
| `eggnog_db` | eggNOG-mapper（5.x，~40+ GB） | **报错**（或 `--skip_eggnog`） |
| `card_db` | RGI/CARD（card.json，GB 级） | **报错**（或 `--skip_rgi`） |
| `pathway_db` | KO→pathway 映射（Phase 14 Pathway 列） | **可选**：列留空 |

read-based 三支缺失时告警跳过（并行旁支，不阻塞主路）；MAG 级阶段缺失时明确
报错（下游依赖其结果，跳过会产生不完整的核心表）—— 设计上**无静默降级、无
伪造产物**。完整清单（体量 / 获取方式 / db_dir 约定树）见
[docs/database.md](database.md)。
