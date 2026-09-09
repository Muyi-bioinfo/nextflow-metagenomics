# 输入文档

## 必需输入

### Samplesheet（CSV）

`--input` 指向一个 CSV 文件，每行一个样本：

```csv
sample,fastq_1,fastq_2,group,batch,host
S01,/data/reads/S01_R1.fastq.gz,/data/reads/S01_R2.fastq.gz,case,batch01,human
S02,/data/reads/S02_R1.fastq.gz,/data/reads/S02_R2.fastq.gz,control,batch01,human
```

| 字段 | 必需 | 说明 |
|------|------|------|
| `sample` | 是 | 样本 ID，同一 samplesheet 内唯一（成为 `meta.id`） |
| `fastq_1` | 是 | R1 FASTQ 路径（绝对路径或相对 launch 目录的相对路径） |
| `fastq_2` | 是* | R2 FASTQ 路径；单端样本留空（`*`：PE 样本必需，SE 样本留空） |
| `group` | 否 | 实验分组（`meta.group`，coassembly V2 预留） |
| `batch` | 否 | 实验批次（`meta.batch`） |
| `host` | 否 | 宿主物种标识（`meta.host`，信息用途） |

### 验证流程（CHECK_SAMPLESHEET process）

1. 校验字段、检查文件存在性，路径相对 launch 目录解析后**补全为绝对路径**；
2. `fastq_2` 为空 → 判定单端（SE），写入 `single_end=True`；否则双端（PE）；
3. 产出 `00_metadata/validated_samplesheet.csv`（含 `single_end` 列），
   主 workflow 据此构建 `tuple(meta, reads)` 通道 —— reads 为 `[R1]`（SE）
   或 `[R1, R2]`（PE）。

FASTQ 文件可为 gzip 压缩或未压缩（由 fastp 处理）。

### 宿主索引（可选）

`--host_index`：宿主基因组的 Bowtie2 索引前缀（不含 `.1.bt2` 后缀，用
`bowtie2-build` 构建）。未提供时**跳过宿主去除**（clean reads 直通下游，
可用 `--skip_host_removal` 显式声明）。

---

## 数据库输入

各分析阶段的数据库经 `--*_db` 参数传入（或经 `--db_dir` 约定树派生），
缺失时的行为（跳过 / 报错）各不相同。完整清单、体量与获取方式见
[database.md](database.md)。

---

## 测试输入（`-profile test`）

`conf/test.config` 定义内置测试运行：

- **合成数据**：`test/data/` 下 2 个 PE 样本（合成宏基因组，总长约 46 kb），
  由 `test/data/make_test_data.py` 生成（可复现）；
- **stub 门卫目录**：`test/data/dummy_dbs/`（read-based 三支的参数守卫为
  `checkIfExists`，stub-run 时指向这些空目录）；
- **pathway 夹具**：`test/data/pathway_test.tsv`（KO→pathway 合成映射，
  含 "ko:" 前缀行）。

> 测试数据总长 46 kb < MetaBAT2 的 200 kb 最小 bin 阈值，test.config 将
> `metabat2_min_bin_size` 降至 10 kb 使 S02 能形成 1 个 bin。测试验证的是
> **全流程可运行与数值交叉一致**，不验证真实宏基因组的生物学结论。
