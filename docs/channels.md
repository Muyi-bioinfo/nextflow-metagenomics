# 通道设计

所有通道使用 `tuple(meta, ...)` 结构，meta 是 Groovy map，贯穿整个工作流
。

## Meta 对象

```groovy
// samplesheet 解析后（Phase 2）
meta = [
    id:         String,   // 样本 ID，同一 samplesheet 内唯一
    group:      String,   // 实验分组
    batch:      String,   // 批次 ID
    host:       String,   // 宿主物种标识
    single_end: Boolean   // false = PE（默认）; true = SE
]

// 组装子工作流中追加（Phase 5）
meta.assembler      = "megahit" | "metaspades"
meta.assembly_mode  = "single" | "coassembly"
meta.samples        = [sample_id, ...]   // single 时长度为 1

// Phase 6/7 继续使用: assembly_id (组装单元 id)
```

SE 样本：reads = `[R1]`（单元素列表）。meta.single_end 控制各模块的命令行构建。

---

## 关键通道一览

| 通道 | 结构 | 来源 | 消费方 |
|------|------|------|--------|
| raw reads | `tuple(meta, [R1, R2])` | samplesheet 解析（Phase 2） | PREPROCESSING |
| clean reads | `tuple(meta, [R1, R2])` | `PREPROCESSING.out.reads` | READ_BASED + ASSEMBLY + MAPPING（三者共用，上游不重复执行） |
| contigs | `tuple(meta, contigs.fa.gz)` | `ASSEMBLY.out.contigs` | MAPPING（Phase 6）、BINNING（Phase 7） |
| depth | `tuple(meta, depth.txt)` | `MAPPING.out.depth` | BINNING（Phase 7） |
| bam | `tuple(meta, bam, bai)` | `MAPPING.out.bam` | ABUNDANCE（Phase 13） |
| MAG records | `tuple(meta, mag_id, mag_fasta)` | `BINNING.out.mags`（Phase 7） | MAG_QC（Phase 8） |
| qualified MAGs | `tuple(meta, mag_id, mag_fasta)` | `MAG_QC.out.qualified_mags`（形状同 mags） | DREP（Phase 9） |
| qc table | `path(mag_qc.tsv)` | `MAG_QC.out.qc_table` | DREP 打分 + INTEGRATION |
| representatives | `tuple(meta, mag_id, mag_fasta)` | `DREP.out.representatives` | TAXONOMY / GENE_PREDICTION / ABUNDANCE / INTEGRATION |
| membership | `path(mag_membership.tsv)` | `DREP.out.membership` | TAXONOMY（成员回填）+ INTEGRATION |
| taxonomy table | `path(mag_taxonomy.tsv)` | `TAXONOMY.out.taxonomy_table` | INTEGRATION + MULTIQC |
| proteins | `tuple(meta, mag_id, proteins.faa)` | `GENE_PREDICTION.out.proteins` | ANNOTATION（Phase 12 三支的唯一输入） |
| genes / gff | `tuple(meta, mag_id, file)` | `GENE_PREDICTION.out.genes / gff` | （发布产物） |
| annotation tables | `path(diamond_hits.tsv)` 等 | `ANNOTATION.out.diamond_table / eggnog_table / rgi_table` | INTEGRATION（键 = (meta_id, mag_id, gene)） |
| abundance | `path(mag_abundance.tsv)` | `ABUNDANCE.out.abundance` | INTEGRATION + MULTIQC |
| integrated | `path(...)` | `INTEGRATION.out.metadata / functional / ...` | MULTIQC（混入 ch_multiqc_files） |
| read-based merged matrices | `path(merged_<level>.tsv)`（glob）/ `path(beta_diversity.tsv)`（optional）/ `path(merged_pathabundance.tsv)` | `READ_BASED_MERGE.out.bracken_merged / beta_diversity / pathabundance_merged` | 跨样本比较（Phase 21 PLOTTING：PCoA / 热图，只消费现成表） |
| multiqc files | 各 Phase QC 产物的混入通道 | Phase 3/4/5/6/8/10/13/14 | MULTIQC（`.collect()` 后传入） |

MAG 级通道的 `mag_id` 格式：`<unit>.<assembler>.<binner>.<bin_no>`（三位补零），
稳定可重现；含组装单元 id（= meta.id，single 模式）。

---

## 数据流（全流程）

```
CHECK_SAMPLESHEET
    ↓ validated_samplesheet.csv
ch_raw_reads:  tuple(meta, [R1, R2])
    ↓
PREPROCESSING
    ├── FASTQC(ch_raw_reads)                          → QC 报告
    ├── FASTP(ch_raw_reads)                           → clean reads + JSON/HTML
    └── HOST_REMOVAL(FASTP.out.reads, ch_host_index)  → nonhost reads
    ↓
PREPROCESSING.out.reads:  tuple(meta, [R1, R2])   ← 单一来源, 三分支共用
    │
    ├──→ READ_BASED(ch_clean_reads)                   (Phase 4, 与 ASSEMBLY 并行)
    │       ├── KRAKEN2 → BRACKEN (× 层级 scatter)     → 03_taxonomy/
    │       └── HUMANN                                 → 04_function/
    │       └──→ READ_BASED_MERGE(bracken_abundance, pathabundance)   (Phase 20)
    │               → 03_taxonomy/combined/ + 04_function/combined/（跨样本矩阵）
    │
    └──→ ASSEMBLY(ch_clean_reads)                     (Phase 5, 与 READ_BASED 并行)
            ├── MEGAHIT / METASPADES (按 params.assembler)
            └── QUAST → ASSEMBLY_SUMMARY              → 05_assembly/
            ↓
ASSEMBLY.out.contigs:  tuple(meta, contigs.fa.gz)
    ├──→ MAPPING(contigs, clean reads)  (Phase 6)     → depth + bam
    └──→ BINNING(contigs, depth)        (Phase 7)     → tuple(meta, mag_id, mag_fasta)
            ↓
MAG_QC(ch_mags)                          (Phase 8)     → qualified_mags + qc_table
    ↓
DREP(qualified_mags, qc_table)           (Phase 9)     → representatives + membership
    ├──→ TAXONOMY(rep, membership)       (Phase 10)    → mag_taxonomy.tsv
    ├──→ GENE_PREDICTION(rep)            (Phase 11)    → proteins
    │       └──→ ANNOTATION(proteins)    (Phase 12)    → 三张键控注释表
    └──→ ABUNDANCE(MAPPING.out.bam, rep) (Phase 13)    → mag_abundance.tsv
            ↓
INTEGRATION(representatives, membership, qc, taxonomy, abundance,
            diamond, eggnog, rgi, assembly_summary)    (Phase 14)
            → 14_integrated/ 两张核心表 + 引用拷贝
            ↓
PLOTTING(bin_summary, qc, taxonomy, membership, abundance,
         bracken_merged, beta_diversity, pathabundance_merged)  (Phase 21)
            → 各 Phase figures/ 7 张 PNG (只消费现成表)
            ↓
MULTIQC(ch_multiqc_files.collect(), multiqc_config)    (Phase 15)
```

---

## 操作符使用场景

### collect — 跨样本聚合（谨慎：ArrayBag 扁平化）

用于需要所有样本数据才能运行的任务：

```groovy
FASTP_SUMMARY(FASTP.out.json.map { m, j -> j }.collect(), parser)
ASSEMBLY_SUMMARY(QUAST.out.tsv.map { m, t -> t }.collect(), parser)
MULTIQC(ch_multiqc_files.collect())
```

**坑**：`.collect()` 聚合的 `tuple(meta, mag_id, fasta)` 会扁平化为 ArrayBag，
process 内无法按三元组还原。集合级批处理（Phase 9/10/12/13/14）统一采用：

- **manifest 类元数据**：`map{...}.collectFile(name:..., sort: true)` 物化为文件；
- **文件列表**：`map{...}.toSortedList()` 传 `path` 输入（返回普通 List，且
  排序保证哈希确定性，`-resume` 可命中）。

### combine — 笛卡尔积 scatter

Bracken 按层级 scatter（一个样本 × N 个层级 = N 个 task）：

```groovy
ch_kraken2_report.combine(ch_bracken_levels)
// → tuple(meta, report, level)
```

### branch — 条件分流

metaSPAdes SE/PE 分流（metaSPAdes 要求双端 library，SE 样本告警跳过）：

```groovy
ch_units.branch {
    paired: !meta.single_end
    single: meta.single_end
}
```

### combine(by:0) — contigs 与 reads 配对（Phase 6）

Phase 6 需要把每个组装单元的 contigs 与「组装它的那些样本的 clean reads」
配对。两侧 meta 形状不同（contigs 侧多了 assembler/assembly_mode/samples），
**不能**按整个 meta join。键落在样本身份上，用 `combine` 而非 `join` ——
`--assembler both` 时同一样本属于两个组装单元，是多对多关系：

```groovy
ch_units_by_sample = ch_index.flatMap { meta, index ->
    meta.samples.collect { sid -> tuple(sid, meta, index) }
}
ch_reads_by_sample = ch_reads.map { meta, reads -> tuple(meta.id, meta, reads) }
ch_units_by_sample.combine(ch_reads_by_sample, by: 0)
```

### groupTuple + groupKey — 按组装单元聚合 BAM（Phase 6）

`groupKey(key, size)` 带上预期分组大小（= 该单元样本数），分组收齐即下发；
键必须带 assembler，否则 `--assembler both` 时两套 contigs 的 BAM 会混进
同一个矩阵：

```groovy
SAMTOOLS_SORT.out.bam
    .map { meta, bam, bai ->
        tuple( groupKey("${meta.assembly_id}::${meta.assembler}", meta.samples.size()), bam, bai )
    }
    .groupTuple()
```

### flatMap — MAG 展开（Phase 7）

SPLIT_BINS 的 output 是 `tuple(meta, path("mags/*.fa"))`（一个 bin 目录 =
N 个文件），emit 块用 flatMap 展开成每条 MAG 一个 tuple（0 bin 时为空通道，
合法）。

---

## 通道多消费

PREPROCESSING.out.reads 同时流入 READ_BASED、ASSEMBLY 与 MAPPING：

- Nextflow 通道可被多个下游 process 同时订阅；
- 上游 task（fastp / host_removal）只执行一次；
- 各分支使用同一份物理文件，不重复计算。

---

## 空通道初始化模式

跳过的分支 emit 空通道，保持 emit 结构完整：

```groovy
ch_contigs = Channel.empty()

if (!params.skip_assembly) {
    ASSEMBLY(ch_reads)
    ch_contigs = ASSEMBLY.out.contigs
}
```

**坑**：`Channel.empty().collect()/collectFile()` 都不发射 —— 空上游时下游
聚合阶段直接不调度任务（不产出空表），上游 skip 已有告警兜底；`ifEmpty(通道)`
会把通道对象本身当值传给下游（DataflowStream 泄漏）—— ifEmpty 只接受
**具体值**，可选表兜底统一用 0 字节哨兵文件（assets/empty.tsv）。

## 可选表与哨兵（Phase 14）

INTEGRATION 的每个输入表都可能因对应 `--skip_*` 为空通道，子工作流内以
0 字节哨兵 ifEmpty 兜底，缺失表的对应列留空、行集合由存在表的并集决定。
哨兵喂给多个输入时需 `stageAs: '<固定名>.tsv'` 规避 Nextflow 输入文件名冲突
（同一哨兵按原名暂存多个输入会报 name collision）。

## Read-based 跨样本合并（Phase 20）

READ_BASED_MERGE 只消费 READ_BASED 的逐样本 emit，产出跨样本宽表矩阵。
输入/输出契约：

```
take:
  ch_bracken_abundance   tuple(meta, level, path)   READ_BASED.out.bracken_abundance
  ch_pathabundance       tuple(meta, path)          READ_BASED.out.pathabundance

emit:
  bracken_merged        path(merged_<level>.tsv)    glob 匹配各层级矩阵（单次任务
                                                   产出全部层级，emit 为文件列表，
                                                   可用 flatten 展开）
  beta_diversity        path(beta_diversity.tsv)    optional（单样本/0 taxa 不发射）
  pathabundance_merged  path(merged_pathabundance.tsv)
```

矩阵 schema（详见 docs/output.md）：

| 矩阵 | 行键 | 列 | 数值 | 缺失 |
|------|------|----|------|------|
| `merged_<level>.tsv` | Bracken `name`（分类单元名） | 样本 ID | `fraction_total_reads`（0-1 相对丰度，原样保留） | 补 `0` |
| `beta_diversity.tsv` | 样本 ID（首格为空，标准距离矩阵格式） | 样本 ID | Bray-Curtis 距离（对称，对角 0） | — |
| `merged_pathabundance.tsv` | HUMAnN pathway 字符串 | 样本 ID | HUMAnN Abundance（RPK，原样保留） | 补 `0` |

聚合模式沿用 ArrayBag 规避：manifest 经 collectFile
物化（sort: true）、文件列表经 toSortedList 传 `path` 输入（声明为 path 输入
保证 -resume 依赖追踪），脚本经 manifest 的 resolve() 双路径打开文件（Phase
12/14 模式）。skip_kraken2 / skip_bracken / skip_humann 各自跳过时对应逐样本
emit 为空通道 → collectFile 不发射 → 对应 merge process 不调度（告警已由
read_based.nf 发出）。合并值取逐样本表的原始字符串（不做重归一化），合并
数值正确性由 bin/merge_read_based.py 单测覆盖（合成表）。

## 结果可视化（Phase 21）

PLOTTING 只消费各 Phase 已产出的 TSV（含 Phase 20 combined 矩阵），产出 7 张
PNG 到对应 Phase 的 `figures/` 子目录（有图才建）。输入/输出契约：

```
take:
  ch_bin_summary          path(bin_summary.tsv)       BINNING.out.summary（漏斗 raw bins 锚）
  ch_qc_table             path(mag_qc.tsv)            MAG_QC.out.qc_table（散点 + 漏斗）
  ch_taxonomy_table       path(mag_taxonomy.tsv)      TAXONOMY.out.taxonomy_table（组成 + 漏斗）
  ch_membership           path(mag_membership.tsv)    DREP.out.membership（漏斗 after dRep）
  ch_abundance            path(mag_abundance.tsv)     ABUNDANCE.out.abundance（丰度热图）
  ch_bracken_merged       path(merged_<level>.tsv)    READ_BASED_MERGE.out.bracken_merged（top taxa）
  ch_beta_diversity       path(beta_diversity.tsv)    READ_BASED_MERGE.out.beta_diversity（PCoA, optional）
  ch_pathabundance_merged path(merged_pathabundance.tsv) READ_BASED_MERGE.out.pathabundance_merged（pathway 热图）

emit:
  versions                path(versions.yml)           仅版本（figures 经 publishDir 发布, 不下游消费）
```

调度语义：每个 plot process 以其输入通道是否为空独立调度 —— 对应 Phase 被 skip
时（通道空）该图不调度、`figures/` 不建。漏斗以 bin_summary 为锚（非空才调度），
其余 3 张可选表（qc/membership/taxonomy）以 0 字节哨兵 `ifEmpty` 兜底（stageAs
固定暂存名规避同名冲突，同 INTEGRATE_METADATA 模式），脚本按「0 字节 = 层级缺失」
跳过该级画剩余漏斗。注意：**不要对 glob emit 的通道用 toSortedList()** —— 空通道
上 toSortedList 会发射空列表（而非不发射），导致 process 以空输入误调度；直接传
通道即可（空通道不调度，单/多文件由 process `path` 输入 + 脚本 `nargs='+'` 接住）。
脚本按列名定位（缺列 SystemExit），空表/单点告警跳过不写 PNG（PNG 输出为
optional），matplotlib Agg backend + 固定 figsize 8×6 / dpi 100（输出 800×600）。

