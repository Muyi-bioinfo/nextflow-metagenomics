// ============================================================================
// MAPPING — reads 比对与覆盖度子工作流 (Phase 6)
//
// ─── 数据流 ───────────────────────────────────────────────────────────────
//
//   contigs ──→ BOWTIE2_BUILD ──→ index ──┐
//                                          ├─→ BOWTIE2_MAP ─→ SAMTOOLS_SORT ─┐
//   clean reads ────────────────────────── ┘   (每样本 × 每组装单元)          │
//                                                                             │
//                        depth.txt ←── CONTIG_DEPTH ←── 按组装单元聚合 BAM ←──┘
//                        (MetaBAT2 --abdFile 直接可用)
//
// ─── 本子工作流的核心是「配对」, 不是比对 ────────────────────────────────
//
// 两个上游通道的 meta 形状**不同**:
//   clean reads   meta = [id, group, batch, host, single_end]
//   contigs       meta = 上面这些 + [assembler, assembly_mode, samples]
//
// 因此 ch_contigs.join(ch_reads) 直接按整个 meta 做键会一条都匹配不上 ——
// 两侧的 map 永远不相等。配对必须落在**样本身份**上:
//
//   组装单元 --(meta.samples)--> 样本 id --(meta.id)--> clean reads
//
// meta.samples 是 assembly.nf 专为此预留的字段 (single 模式下 = [meta.id])。
// 用它而不是直接用 meta.id 配对, 是为了 coassembly 接入时本文件无需改动:
// 届时一个组装单元的 samples 有多个, flatMap 自然展开成多个比对任务。
//
// 用 combine(by:0) 而非 join(by:0): --assembler both 时同一个样本 id 会同时
// 出现在 megahit 与 metaspades 两个组装单元中, 是**多对多**关系; join 要求
// 键在两侧唯一, 会丢数据。
//
// ─── 比对任务的 meta ──────────────────────────────────────────────────────
//
// 一个比对任务同时属于「一个样本」和「一套 contigs」, 两者的身份都必须保留:
//   meta.id           比对任务 id, 用于文件名 (如 S01.megahit)
//   meta.sample       reads 来自哪个样本
//   meta.assembly_id  比对到哪个组装单元
//   meta.assembler / assembly_mode / samples   组装单元的其余身份信息
// 其余字段 (group/batch/host/single_end) 取**样本**的值 —— 这条记录描述的是
// 「这个样本的 reads」, 而组装单元的 group 在 coassembly 下可能是聚合值。
//
// ─── 输出目录 ─────────────────────────────────────────────────────────────
//   06_mapping/logs/    bowtie2 日志 + flagstat/idxstats/stats
//   06_mapping/depth/   <unit>.<assembler>.depth.txt   ← Phase 7 MetaBAT2 输入
//   06_mapping/bam/     排序 BAM + BAI            (--save_bam)
//   06_mapping/index/   Bowtie2 索引              (--save_bowtie2_index)
// ============================================================================

include { BOWTIE2_BUILD } from '../../modules/local/mapping/bowtie2_build.nf'
include { BOWTIE2_MAP   } from '../../modules/local/mapping/bowtie2_map.nf'
include { SAMTOOLS_SORT } from '../../modules/local/mapping/samtools_sort.nf'
include { CONTIG_DEPTH  } from '../../modules/local/mapping/coverage.nf'

// ---------------------------------------------------------------------------
// 组装单元的唯一键。
// 必须带上 assembler: --assembler both 时同一样本有两套独立 contigs, 只用
// id 做键会把两者的 BAM 混进同一个覆盖度矩阵。
// ---------------------------------------------------------------------------
def unit_key(String id, String assembler) {
    return "${id}::${assembler}".toString()
}

workflow MAPPING {

    take:
    ch_contigs   // channel: [ val(meta), path(contigs) ]  meta 含 assembler/assembly_mode/samples
    ch_reads     // channel: [ val(meta), path(reads) ]    reads=[R1,R2] (PE) 或 [R1] (SE)

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // =====================================================================
    // 1. 每个组装单元建一次 Bowtie2 索引
    //    参考序列来自通道 (本次组装的 contigs), 没有任何外部路径参与
    // =====================================================================
    BOWTIE2_BUILD(ch_contigs)
    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions.first())

    // =====================================================================
    // 2. 配对 (组装单元索引) × (该单元包含的每个样本的 clean reads)
    // =====================================================================
    ch_reads_by_sample = ch_reads
        .map { meta, reads -> tuple(meta.id, meta, reads) }

    ch_units_by_sample = BOWTIE2_BUILD.out.index
        .flatMap { meta, index ->
            // single: samples=[id] → 展开成 1 条; coassembly: 展开成 N 条
            meta.samples.collect { sample_id -> tuple(sample_id, meta, index) }
        }

    ch_map_input = ch_units_by_sample
        .combine(ch_reads_by_sample, by: 0)
        .map { sample_id, asm_meta, index, read_meta, reads ->
            // single 模式下组装单元就是样本本身, 文件名无需重复样本 id;
            // coassembly 下单元 id 与样本 id 不同, 写成 <sample>.vs.<unit>.<asm>
            def pair_id = (read_meta.id == asm_meta.id)
                ? "${read_meta.id}.${asm_meta.assembler}"
                : "${read_meta.id}.vs.${asm_meta.id}.${asm_meta.assembler}"

            // 以 asm_meta 为底再覆盖样本字段: 组装单元的字段一个都不丢
            def map_meta = asm_meta + [
                id:          pair_id.toString(),
                assembly_id: asm_meta.id,
                sample:      read_meta.id,
                group:       read_meta.group,
                batch:       read_meta.batch,
                host:        read_meta.host,
                single_end:  read_meta.single_end
            ]

            tuple(map_meta, reads, index)
        }

    BOWTIE2_MAP(ch_map_input)
    ch_versions = ch_versions.mix(BOWTIE2_MAP.out.versions.first())

    // bowtie2 日志含总体比对率, MultiQC 的 bowtie2 模块可直接解析
    ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_MAP.out.log.map { meta, log -> log })

    // =====================================================================
    // 3. 排序 + 索引 + 比对统计
    //    覆盖度计算要求坐标排序 BAM, 这一步是硬前提而非可选项
    // =====================================================================
    SAMTOOLS_SORT(BOWTIE2_MAP.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    ch_multiqc_files = ch_multiqc_files
        .mix(SAMTOOLS_SORT.out.flagstat.map { meta, f -> f })
        .mix(SAMTOOLS_SORT.out.idxstats.map { meta, f -> f })
        .mix(SAMTOOLS_SORT.out.stats.map    { meta, f -> f })

    // =====================================================================
    // 4. 按组装单元聚合 BAM → contig 覆盖度矩阵
    //
    // 聚合的是「比对到同一套 contigs 的全部样本 BAM」—— 深度矩阵的列就是
    // 这些样本。single 模式下每单元 1 个 BAM; coassembly 下为该组全部样本,
    // MetaBAT2 才能用上跨样本共变异信号。
    //
    // groupKey(key, size) 携带预期分组大小 (= 该单元的样本数), 使 groupTuple
    // 在收齐即可下发, 不必等整条通道结束。
    // =====================================================================
    ch_bams_by_unit = SAMTOOLS_SORT.out.bam
        .map { meta, bam, bai ->
            def key = unit_key(meta.assembly_id, meta.assembler)
            tuple( groupKey(key, meta.samples.size()), bam, bai )
        }
        .groupTuple()
        .map { key, bams, bais ->
            // 按文件名排序: groupTuple 的到达顺序不确定, 而 BAM 顺序既决定
            // 深度矩阵的列顺序, 也进入 task hash —— 不排序会破坏 -resume
            tuple( key.toString(), bams.sort { it.name }, bais.sort { it.name } )
        }

    // 组装单元的原始 meta 从 contigs 通道取回, 不做任何重建
    ch_units = ch_contigs
        .map { meta, contigs -> tuple( unit_key(meta.id, meta.assembler), meta ) }

    ch_depth_input = ch_units
        .join(ch_bams_by_unit)
        .map { key, meta, bams, bais -> tuple(meta, bams, bais) }

    CONTIG_DEPTH(ch_depth_input)
    ch_versions = ch_versions.mix(CONTIG_DEPTH.out.versions.first())

    emit:
    // ---- 比对结果 ----
    bam           = SAMTOOLS_SORT.out.bam        // [ val(meta), path(bam), path(bai) ]
    flagstat      = SAMTOOLS_SORT.out.flagstat   // [ val(meta), path(flagstat) ]
    idxstats      = SAMTOOLS_SORT.out.idxstats   // [ val(meta), path(idxstats) ]

    // ---- 覆盖度 ----
    // [ val(meta), path(depth.txt) ]  meta = 组装单元 meta (与 contigs 通道一致)
    // Phase 7: METABAT2 以 unit_key 与 contigs 通道 join 即可配对
    depth         = CONTIG_DEPTH.out.depth

    // ---- 汇总 ----
    multiqc_files = ch_multiqc_files             // path(*)  bowtie2 日志 + samtools 统计
    versions      = ch_versions                  // path(versions.yml)
}
