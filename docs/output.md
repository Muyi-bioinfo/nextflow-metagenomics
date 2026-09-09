# 输出文档

每次 workflow run 的输出写到 `results/<batch_id>/`（batch 目录），内部使用固定
编号目录。目录树与各 process 的 `publishDir` 一一对应 —— 被 `--skip_*` 跳过的
阶段不产生对应目录，`--save_*` 关闭的中间产物不发布（默认可再生的不发布）。

```text
results/<batch_id>/
├── 00_metadata/
├── 01_qc/
├── 02_host_removal/
├── 03_taxonomy/
├── 04_function/
├── 05_assembly/
├── 06_mapping/
├── 07_binning/
├── 08_mag_qc/
├── 09_dereplication/
├── 10_mag_taxonomy/
├── 11_gene_prediction/
├── 12_annotation/
├── 13_abundance/
├── 14_integrated/
├── 99_multiqc/
└── run_reports/              (可选, -with-* 显式指定路径)
```

---

## 00_metadata/ — 输入验证

| 文件 | 说明 |
|------|------|
| `validated_samplesheet.csv` | 验证后的 samplesheet（补全绝对路径/默认值，列含 `single_end`） |

## 01_qc/ — 预处理质控

| 文件 | 说明 |
|------|------|
| `fastqc/<sample>_fastqc.{html,zip}` | FastQC 报告（raw reads） |
| `fastp/<sample>.fastp.{json,html}` | fastp 报告 |
| `fastp/clean_reads/<sample>.fastp.fastq.gz` | fastp 清洁 reads（`--save_trimmed`，默认 true） |
| `fastp_summary.tsv` | 逐样本 fastp 汇总（解析自 JSON，MultiQC 输入） |

## 02_host_removal/ — 宿主去除

| 文件 | 说明 |
|------|------|
| `logs/<sample>.{log,flagstat,idxstats}` | Bowtie2 宿主比对日志 + samtools 统计 |
| `reads/<sample>_nonhost_R*.fastq.gz` | 去宿主后 reads（`--save_host_removed`，默认 false，体积大） |

未提供 `--host_index` 时跳过宿主去除，clean reads 直通下游。

## 03_taxonomy/ — Read-based 物种分类（Kraken2 → Bracken）

| 文件 | 说明 |
|------|------|
| `kraken2/<sample>.kraken2.report.txt` | Kraken2 标准报告（Bracken 输入，MultiQC 可解析） |
| `kraken2/<sample>.kraken2.log` | Kraken2 运行日志 |
| `kraken2/classifications/<sample>.kraken2.output.txt.gz` | 逐 read 分类结果（`--save_kraken2_output`，默认 false，可达数 GB） |
| `bracken/<sample>.bracken.<level>.tsv` | Bracken 丰度表（每层级一个，层级见 `--bracken_levels`） |
| `bracken/<sample>.bracken.<level>.report.txt` | Bracken 重估后的 Kraken 风格报告 |

## 04_function/ — Read-based 功能谱（HUMAnN）

| 文件 | 说明 |
|------|------|
| `humann/<sample>_genefamilies.tsv` | 基因家族丰度表 |
| `humann/<sample>_pathabundance.tsv` | 通路丰度表 |
| `humann/<sample>_pathcoverage.tsv` | 通路覆盖度表 |
| `humann/<sample>.humann.log` | HUMAnN 运行日志 |

## 05_assembly/ — 组装与质控

| 文件 | 说明 |
|------|------|
| `megahit/<unit>/final.contigs.fa.gz` 等 | MEGAHIT 输出目录（`<unit>` = 组装单元 ID） |
| `metaspades/<unit>/` | metaSPAdes 输出目录（`--assembler metaspades|both`） |
| `quast/<unit>/` | QUAST 报告目录（含 `report.tsv`、`report.html`） |
| `assembly_summary.tsv` | 逐组装 QUAST 汇总表（`--skip_quast` 时无此文件） |

## 06_mapping/ — Reads 比对与覆盖度

| 文件 | 说明 |
|------|------|
| `logs/` | Bowtie2 比对日志 + samtools flagstat/idxstats/stats |
| `depth/<unit>.<assembler>.depth.txt` | MetaBAT2 兼容覆盖度矩阵（Phase 7 输入） |
| `bam/<sample>.sorted.bam{,.bai}` | 排序 BAM + 索引（`--save_bam`，默认 false，可达数十 GB） |
| `index/<unit>.<assembler>/` | contigs 的 Bowtie2 索引（`--save_bowtie2_index`，默认 false） |

## 07_binning/ — MAG 分箱

| 文件 | 说明 |
|------|------|
| `metabat2/<unit>.<assembler>/` | MetaBAT2 原始 bins 目录 |
| `logs/` | MetaBAT2 日志 |
| `mags/<unit>.<assembler>.<binner>.<bin_no>.fa` | 提取出的 MAG FASTA（MAG ID 稳定可重现） |
| `<unit>.<assembler>.<binner>_summary.tsv` | 逐组装单元 bin 统计表 |
| `bin_summary.tsv` | 全局 bin 汇总（列：mag_id/assembly_unit/assembler/assembly_mode/binner/source_bin/n_contigs/total_bp/...） |

## 08_mag_qc/ — MAG 质控（CheckM2）

| 文件 | 说明 |
|------|------|
| `checkm2/<mag_id>.qc.tsv` | 逐 MAG CheckM2 结果 |
| `mag_qc.tsv` | 全局 QC 汇总表（完整度/污染度，dRep 打分输入） |

## 09_dereplication/ — MAG 去冗余（dRep）

| 文件 | 说明 |
|------|------|
| `dereplicated_genomes/` | 代表 MAG FASTA 目录 |
| `data_tables/` | dRep 聚类原始数据表（Cdb.csv / Wdb.csv 等） |
| `drep_log/` | dRep 日志 |
| `mag_membership.tsv` | 成员表：`mag_id \t sample \t representative_mag_id`（无表头） |

## 10_mag_taxonomy/ — MAG 分类（GTDB-Tk）

| 文件 | 说明 |
|------|------|
| `gtdbtk/` | 平铺的 summary TSV（bac120/ar53） |
| `gtdbtk_out/` | GTDB-Tk 原始输出目录 |
| `mag_taxonomy.tsv` | 分类汇总表（sample/mag_id/rep_mag_id/domain..species，成员 MAG 分类由代表回填） |

## 11_gene_prediction/ — 基因预测（Prodigal）

| 文件 | 说明 |
|------|------|
| `<mag_id>.genes.fna` | 预测基因核酸序列 |
| `<mag_id>.proteins.faa` | 预测蛋白序列（Phase 12 输入） |
| `<mag_id>.gff` | GFF3 注释 |

## 12_annotation/ — 功能注释（三支并行）

| 文件 | 说明 |
|------|------|
| `diamond/<mag_id>.diamond.tsv` | 逐 MAG DIAMOND blastp 原始命中 |
| `diamond/diamond_hits.tsv` | 汇总表（meta_id/mag_id/gene/sseqid/pident/length/evalue/bitscore，每 gene 最优命中） |
| `eggnog/<mag_id>.emapper.{annotations,seed_orthologs}` | 逐 MAG eggNOG-mapper 原始注释 |
| `eggnog/eggnog_annotations.tsv` | 汇总表（meta_id/mag_id/gene/COG_category/Description/Preferred_name/EC/KO/GO/KEGG_Pathway） |
| `rgi/<mag_id>.rgi.json` | 逐 MAG RGI 原始结果 |
| `rgi/rgi_annotations.tsv` | 汇总表（meta_id/mag_id/gene/ARO/ARO_accession/AMR_gene_family/drug_class/resistance_mechanism/pct_identity/model_type） |

## 13_abundance/ — MAG 丰度（CoverM）

| 文件 | 说明 |
|------|------|
| `mag_abundance.tsv` | 跨样本丰度矩阵：行 = mag_id、列 = 样本 ID，相对丰度 0-1 |

## 14_integrated/ — 整合结果（核心交付）

| 文件 | 说明 |
|------|------|
| `mag_metadata.tsv` | **核心表**：MAG_ID / Sample / Completeness / Contamination / Genome_size / GC / GTDB_taxonomy / Species + 各样本丰度列（对应表被 skip 时对应列留空，未检出样本补 0） |
| `mag_functional_annotation.tsv` | **核心表**：MAG_ID / Gene / KO / COG / GO / Pathway / ARG（三表 gene 并集；Pathway 需可选 `--pathway_db`，缺失时留空） |
| `mag_membership.tsv` | 成员表引用拷贝（内容同 09_dereplication） |
| `assembly_summary.tsv` | 组装汇总引用拷贝（`--skip_quast`/`--skip_assembly` 时无此文件） |

## 99_multiqc/ — MultiQC 汇总报告

| 文件 | 内容 |
|------|------|
| `multiqc_report.html` | 全 batch 汇总报告 |
| `multiqc_report_data/` | 报告底层数据（各模块解析结果） |
| `versions.yml` | 各工具版本 |

覆盖模块：fastqc / fastp / bowtie2（宿主去除 + 比对）/ quast，以及自定义 MAG 级
表 mag_qc.tsv / mag_taxonomy.tsv / mag_metadata.tsv / mag_functional_annotation.tsv /
mag_abundance.tsv（MultiQC custom content 表格式解析，配置见
`modules/local/qc/multiqc_config.yaml`）。

对应 Phase 被 skip 或数据库不可用时，缺失的模块/表不出现于报告（不伪造）。
跳过 `--skip_multiqc` 则不生成整个 99_multiqc/。

## run_reports/ — Nextflow 运行报告（可选）

`-with-*` 为 Nextflow 内建参数，产物路径显式指定到
`results/<batch_id>/run_reports/`：

```bash
nextflow run main.nf -profile test \
    -with-report   results/<batch_id>/run_reports/report.html \
    -with-timeline results/<batch_id>/run_reports/timeline.html \
    -with-trace    results/<batch_id>/run_reports/trace.txt \
    -with-dag      results/<batch_id>/run_reports/dag.svg
```

| 参数 | 产物 | 说明 |
|------|------|------|
| `-with-report` | report.html | 运行资源使用与任务统计 |
| `-with-timeline` | timeline.html | 任务时间线（甘特图） |
| `-with-trace` | trace.txt | 逐任务 trace 记录 |
| `-with-dag` | dag.svg | 流程 DAG（需系统安装 graphviz） |

---

## 发布规则

- **batch 目录 ≠ work 目录**：`results/<batch_id>/` 存放业务结果；`work/`
  存放 Nextflow task cache，`-resume` 依赖 work/ 而非 batch 目录。
- **默认不发布的中间产物**（`--save_*` 开关控制）：fastp 清洁 reads
  （`--save_trimmed`，默认开）、去宿主 reads（`--save_host_removed`）、
  排序 BAM（`--save_bam`）、Bowtie2 索引（`--save_bowtie2_index`）、
  Kraken2 逐 read 输出（`--save_kraken2_output`）、metaSPAdes assembly graph
  （`--save_assembly_graph`）、未分箱序列（`--save_unbinned`）。
