# 数据库文档

本文档描述 nextflow-metagenomics 使用的所有数据库与参考文件：是什么、大致体量、
经哪个参数传入、被哪个 Phase 使用、**缺失时的行为**（报错 / 跳过 / 留空）。
所有路径均通过 `params.*` 传入，仓库中不出现任何硬编码路径。

---

## 数据库清单

| 数据库 | 用途 | Phase | 参数 | 大致体量 | 缺失时行为 |
|--------|------|-------|------|----------|-----------|
| 宿主基因组索引 | 宿主去除 (Bowtie2) | 3 | `--host_index` | 视宿主基因组而定（人类约 3-4 GB） | **跳过宿主去除**：clean reads 直通下游（可显式 `--skip_host_removal`） |
| Kraken2 标准库 | read 级物种分类 | 4 | `--kraken2_db` | 标准库磁盘数十 GB 起；`--kraken2_memory_mapping false`（默认）时整库入内存，约 40-70 GB RAM | **告警并跳过 Kraken2/Bracken 分支**（不报错，不虚构） |
| Bracken kmer_distrib | read 丰度估计 | 4 | `--bracken_db` | 内嵌于 Kraken2 库内（无独立体量） | 默认 `bracken_db ?: kraken2_db` 自动跟随 Kraken2 库；两者皆缺时随 Kraken2 分支告警跳过 |
| HUMAnN 三库 | 功能谱分析 | 4 | `--humann_db` / `--humann_nucleotide_db` / `--humann_protein_db` / `--metaphlan_db` | ChocoPhlAn ~5 GB、UniRef ~6-14 GB（uniref50/uniref90）、MetaPhlAn ~1-2 GB | **告警并跳过 HUMAnN 分支**（不报错，不虚构） |
| CheckM2 DIAMOND DB | MAG 完整度/污染度 | 8 | `--checkm2_db` | 约 3 GB | **明确报错**（不静默跳过）；或 `--skip_mag_qc` |
| GTDB-Tk 参考库 | MAG 物种分类 | 10 | `--gtdbtk_db` | R220+ 解压后约 110 GB | **明确报错**；或 `--skip_taxonomy` |
| dRep (MASH/FastANI) | MAG 去冗余 | 9 | 无（无 `drep_db` 参数） | 随 dRep 环境自带 | 无外部依赖；质量打分复用 Phase 8 的 QC 表（`--skip_mag_qc` 必须连带 `--skip_dereplication`，由 DREP 子工作流守卫报错） |
| NR DIAMOND 索引 | 蛋白比对注释 | 12 | `--diamond_db` | NR `.dmnd` 数十 GB；`makedb` 需 100+ GB 内存 | **明确报错**；或 `--skip_diamond` / `--skip_annotation` |
| eggNOG 数据目录 | 功能注释 | 12 | `--eggnog_db` | eggNOG 5.x 约 40+ GB（`eggnog.db` + `eggnog_proteins.dmnd`） | **明确报错**；或 `--skip_eggnog` / `--skip_annotation` |
| CARD card.json | 抗性基因注释 | 12 | `--card_db` | GB 级 | **明确报错**；或 `--skip_rgi` / `--skip_annotation` |
| Prodigal | 基因预测 | 11 | 无 | 独立二进制 | 无外部依赖 |
| CoverM | MAG 丰度 | 13 | 无 | 自带比对器 | 无外部依赖（复用 Phase 6 排序 BAM） |
| KO→pathway 映射 | 功能表 Pathway 列 | 14 | `--pathway_db` | 文本级（MB 级） | **可选**：未提供 → Pathway 列留空（不虚构）；提供了路径但文件不存在 → checkIfExists 明确报错 |

### 缺失行为小结

- **读段级分支（Kraken2/Bracken/HUMAnN，Phase 4）**：数据库缺失时**告警跳过**，
  不影响其他分支 —— 这是设计选择：read-based 是并行旁支，组装/MAG 主路不依赖它。
- **MAG 级阶段（CheckM2/GTDB-Tk/DIAMOND/eggNOG/RGI，Phase 8/10/12）**：数据库
  缺失时**明确报错**，不静默降级 —— 下游（dRep 打分、成员回填、整合 join）依赖
  这些结果，跳过会产生伪造或不完整的核心表。要跳过需显式加对应 `--skip_*`。
- **`--pathway_db`（Phase 14）**：唯一「可选」外部依赖 —— 未提供是合法空值
  （Pathway 列留空），路径错拼是错误（报错）。两者不可混淆。
- **宿主索引（Phase 3）**：未提供时跳过宿主去除，clean reads 直通 —— 测试数据
  与无宿主样本的正常用法。

---

## 标准目录布局（Phase 17：`params.db_dir`）

`db_dir` 非空时，各工具数据库按约定树**派生默认值**，一条参数替代八条：

```
<db_dir>/
├── kraken2/                  # Kraken2 标准库（hash.k2d 等）
│                             #   + bracken 的 kmer_distrib（随库构建产出）
├── humann/                   # HUMAnN 父目录（humann_databases 布局）
│   ├── chocophlan/           #   ChocoPhlAn 核酸数据库
│   ├── uniref/               #   UniRef 蛋白数据库（含 .dmnd）
│   └── metaphlan/            #   MetaPhlAn 数据库（prescreen 用，含 *.pkl）
├── checkm2/                  # CheckM2 DIAMOND 数据库
├── gtdbtk/                   # GTDB-Tk 参考库（release 目录）
├── diamond/                  # NR 蛋白库 .dmnd
├── eggnog/                   # eggNOG 数据目录（eggnog.db 等）
├── card/                     # CARD card.json
└── pathway/                  # KO→pathway 映射文件（**文件节点**，非目录）
```

### 派生与优先级规则

1. **`db_dir` 本身**：命令行 `--db_dir` > gitignored 的 `conf/local.config`
   （存在时自动读入）> 默认 `null`。
2. **逐工具派生**：`<db_dir>/<工具>` 子目录**存在**时取之；子目录缺失或
   `db_dir` 为 null 时该参数维持 `null` —— 现有守卫照常报错，不静默跳过。
3. **显式优先**：命令行/profile 显式传入的 `--gtdbtk_db /custom/path` 等
   覆盖派生值（配置求值在 CLI 参数应用之后）。
4. **bracken 复用 kraken2 库**：`--bracken_db` 不在派生清单中 ——
   read_based.nf 内 `bracken_db ?: kraken2_db` 使其自动跟随 Kraken2 库
   （kmer_distrib 通常就建在 Kraken2 库目录内），显式 `--bracken_db` 仍优先。
5. **chocophlan/uniref 跟随 humann 父目录**：`--humann_nucleotide_db` /
   `--humann_protein_db` 按 `<humann_db>/chocophlan`、`<humann_db>/uniref`
   推导（子目录存在时采用），显式传参逐个覆盖；MetaPhlAn 库必须显式给出
   （或经 `<db_dir>/humann/metaphlan` 派生）。
6. **配置求值顺序实测**（Nextflow 26.04）：配置解析器禁止 if/try 语句，
   上述逻辑以 ConfigSlurper + 顶层赋值表达式实现；CLI 参数在配置求值时已
   可见，故优先级语义成立。四场景实测（派生全链路 / 缺子目录报错 / 显式
   覆盖 / local.config 自动读入）均已通过。

### 本地约定目录（示例，非仓库内容）

本地以 `/data/databases/MAG-db/` 为 `db_dir` 根（与 mNGS-db /
tNGS-db 同根约定），经 gitignored 的 `conf/local.config` 提供 —— 仓库中
不出现个人路径，CI/他人环境无此文件时 `db_dir` 维持 null，行为不变。

---

## 各数据库获取方式

| 数据库 | 获取方式 |
|--------|----------|
| 宿主基因组索引 | `bowtie2-build <host_genome.fa> <prefix>`，`--host_index` 传前缀（不含 `.1.bt2` 后缀） |
| Kraken2 标准库 | `kraken2-build --standard --db <dir>`（或官方预构建 tarball） |
| Bracken kmer_distrib | 对 Kraken2 库运行 `bracken-build`（需知道构建时读长，`--bracken_read_length` 必须与之匹配） |
| HUMAnN 三库 | `humann_databases --download chocophlan/uniref/utility_mapping <dir>` + MetaPhlAn 数据库（`metaphlan --install`） |
| CheckM2 DB | `checkm2 database --download --path <dir>`（~3 GB DIAMOND DB） |
| GTDB-Tk 参考库 | GTDB-Tk 2.x 兼容版本（R220+），官方 `download-db.sh` 或 tarball，解压后约 110 GB |
| NR DIAMOND 索引 | NCBI NR 蛋白库 `diamond makedb --in nr.faa -d nr.dmnd`（makedb 需 100+ GB 内存） |
| eggNOG 数据目录 | eggNOG-mapper 官方 `download_eggnog_data.py`（`eggnog.db` + `eggnog_proteins.dmnd`） |
| CARD card.json | CARD 官网下载；RGI_LOAD 一次性 `rgi load --card_json <card.json> --local` 载入本地库 |
| KO→pathway 映射 | 两列 TSV（`KO \t pathway`，无表头；KO 列可带 `ko:` 前缀，加载时剥离；同一 KO 多行 = 多 pathway），如 KEGG 的 ko_pathway 映射数据 |

## 模块与容器的数据库访问

- **conda 环境（nf-meta）**：数据库路径按 `--*_db` 传入即可。
- **容器模式（docker/singularity）**：路径参数会随 task 挂载；RGI 的
  `rgi main` 不接收数据库路径参数 —— 卡库由 RGI_LOAD process 一次性
  `rgi load --card_json ... --local` 载入本地库（conda 环境内直接可用；
  容器部署时可在镜像内预载入，`--card_db` 仍强制提供作守卫）。
