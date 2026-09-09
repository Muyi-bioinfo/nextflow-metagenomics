// ============================================================================
// BINNING — MAG 分箱子工作流 (Phase 7)
//
// ─── 数据流 ───────────────────────────────────────────────────────────────
//
//   contigs ────┐
//               ├──→ join(unit_key) ──→ METABAT2 ──→ bins ──→ SPLIT_BINS ──→ MAGs
//   depth ──────┘                                                  │
//                                                                  └──→ BIN_SUMMARY
//
// ─── 与 Phase 6 (mapping) 的配对 ──────────────────────────────────────────
//
// 两个上游通道:
//   ch_contigs   [ val(meta), path(contigs) ]   meta 含 id/assembler/assembly_mode/samples
//   ch_depth     [ val(meta), path(depth.txt) ] 同上 (与 contigs 的 meta 一致)
//
// 配对键 = unit_key(meta.id, meta.assembler) = "${id}::${assembler}"
// 必须带上 assembler: --assembler both 时同一样本有两套独立 contigs, 只用
// id 做键会把两者的 depth 混给错的 contigs。
//
// ─── SPLIT_BINS 的输出通道展开 ───────────────────────────────────────────
//
// SPLIT_BINS 产出 `tuple(meta, path("mags/*.fa"))`; bins 目录有 N 个 bin 时
// 就是 N 个文件。下游 Phase 8 (CheckM2) 要逐 MAG 处理, 因此 emit 块用
// flatMap 展开成 `tuple(meta, mag_id, mag_fasta)` —— 每个 MAG 一条记录。
//
// MAG ID 从文件名提取: <unit>.<assembler>.<binner>.<NNN>.fa → 去掉 .fa 后缀。
//
// ─── 输出 ─────────────────────────────────────────────────────────────────
//   07_binning/metabat2/<unit>.<asm>/bins/   MetaBAT2 原始 bin 目录
//   07_binning/mags/                         全局唯一命名的 MAG FASTA
//   07_binning/<unit>.<asm>.<binner>_summary.tsv   per-unit 统计表
//   07_binning/bin_summary.tsv               全局汇总表 (Phase 14 整合结果依赖此表)
//   07_binning/logs/                         MetaBAT2 日志
// ============================================================================

include { METABAT2     } from '../../modules/local/binning/metabat2.nf'
include { SPLIT_BINS   } from '../../modules/local/binning/split_bins.nf'
include { BIN_SUMMARY  } from '../../modules/local/binning/bin_summary.nf'

// ---------------------------------------------------------------------------
// 组装单元的唯一键 (与 subworkflows/local/mapping.nf 定义一致)
// ---------------------------------------------------------------------------
def unit_key(String id, String assembler) {
    return "${id}::${assembler}".toString()
}

workflow BINNING {

    take:
    ch_contigs   // channel: [ val(meta), path(contigs) ]  meta 含 id/assembler/assembly_mode/samples
    ch_depth     // channel: [ val(meta), path(depth.txt) ]  同上

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // =====================================================================
    // 1. 配对 contigs + depth (按组装单元)
    // =====================================================================
    ch_contigs_keyed = ch_contigs
        .map { meta, contigs -> tuple(unit_key(meta.id, meta.assembler), meta, contigs) }

    ch_depth_keyed = ch_depth
        .map { meta, depth -> tuple(unit_key(meta.id, meta.assembler), depth) }

    ch_binning_input = ch_contigs_keyed
        .join(ch_depth_keyed)
        .map { key, meta, contigs, depth -> tuple(meta, contigs, depth) }

    // =====================================================================
    // 2. MetaBAT2 分箱
    // =====================================================================
    METABAT2(ch_binning_input)
    ch_versions = ch_versions.mix(METABAT2.out.versions.first())

    // MetaBAT2 日志可被 MultiQC 解析 (如果有对应模块), 暂不确定是否需要
    // ch_multiqc_files = ch_multiqc_files.mix(METABAT2.out.log)

    // =====================================================================
    // 3. 拆分 bin 目录 → 独立 MAG FASTA + per-unit summary
    // =====================================================================
    ch_splitter = file("${projectDir}/bin/split_bins.py", checkIfExists: true)

    SPLIT_BINS(METABAT2.out.bins, ch_splitter)
    ch_versions = ch_versions.mix(SPLIT_BINS.out.versions.first())

    // =====================================================================
    // 4. 展开成 tuple(meta, mag_id, mag_fasta) —— 每个 MAG 一条记录
    //
    // SPLIT_BINS.out.mags = [ val(meta), path("mags/*.fa") ]
    // 0 个 bin 时是 [ meta, [] ]; flatMap 展开后自然是空通道, 不触发下游。
    // =====================================================================
    ch_mags = SPLIT_BINS.out.mags
        .flatMap { meta, mag_files ->
            // mag_files 可能是单个 Path (1 个 MAG) 或 List<Path> (N 个 MAG) 或空列表
            def files = mag_files instanceof List ? mag_files : (mag_files ? [mag_files] : [])
            files.collect { mag_fasta ->
                // MAG ID = 文件名去掉 .fa 后缀
                def mag_id = mag_fasta.name.replaceAll(/\.fa$/, '')
                tuple(meta, mag_id, mag_fasta)
            }
        }

    // =====================================================================
    // 5. 汇总 per-unit summary → bin_summary.tsv
    // =====================================================================
    BIN_SUMMARY(SPLIT_BINS.out.summary.collect())
    ch_versions = ch_versions.mix(BIN_SUMMARY.out.versions)

    emit:
    // ---- MAG 输出 ----
    // [ val(meta), val(mag_id), path(mag.fa) ]  每个 MAG 一条; Phase 8 CheckM2 从此接入
    mags          = ch_mags

    // ---- 中间产物 ----
    bins          = METABAT2.out.bins         // [ val(meta), path(bins_dir) ]  原始 bin 目录 (调试用)
    summary       = BIN_SUMMARY.out.tsv       // path(bin_summary.tsv)  全局汇总表

    // ---- 汇总 ----
    multiqc_files = ch_multiqc_files          // path(*)  暂为空
    versions      = ch_versions               // path(versions.yml)
}
