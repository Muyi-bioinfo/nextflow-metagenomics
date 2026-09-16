// ============================================================================
// ABUNDANCE — MAG 丰度子工作流 (Phase 13)
//
// 输入:  bams             = MAPPING.out.bam  [ val(meta), path(bam), path(bai) ]
//                           Phase 6 坐标排序 + 已建索引的 BAM (clean reads)
//        representatives  = DREP.out.representatives
//                           [ val(meta), path(mag_id), path(mag_fasta) ]
//
// 处理:  集合级批处理 —— 全部 BAM + 全部代表 MAG 聚合成一次 CoverM genome
//        调用, 直接产出跨样本丰度矩阵 (行 = mag_id, 列 = 样本)。
//        BAM 直接复用 Phase 6 产物: MAG FASTA 的 contig header 原样保留
//        (bin/split_bins.py), 与 BAM 的 contig 名一致, 无需重新比对。
//
// 输出:  abundance = path(mag_abundance.tsv)   矩阵键控:
//                                                行 = mag_id, 列 = 样本 ID
//
// ─── 列名规则 ─────────────────────────────────────────────────────────────
// 单组装器 (默认): 列名 = 样本 ID (BAM 与样本一一对应)。
// --assembler both: 同一样本有两套 BAM (各对应一套 contigs), 列名消歧为
// <sample>.<assembler> —— 两套 contigs 的丰度不可合并为一列, 且 MAG 的
// mag_id 本身带 assembler 字段, 消歧列名与其一一对应。
//
// ─── 聚合规则 (ArrayBag 规避, 详见文档) ───────────────────
//   - BAM 清单: map{...}.collectFile(name:..., sort: true) 物化为文件;
//   - BAM/BAI 与 MAG FASTA: map{...}.toSortedList() 传 path 输入
//     (返回普通 List 而非 ArrayBag, 排序保证 -resume 哈希稳定)。
// 空通道 (--skip_mapping / --skip_binning 链式传导) 时 collectFile /
// toSortedList 都不发射, COVERM 不调度, 本阶段正常收尾不报错 —— 上游
// skip 守卫 (mag.nf) 已有明确告警兜底。
// ============================================================================

include { COVERM } from '../../modules/local/abundance/coverm.nf'

workflow ABUNDANCE {

    take:
    ch_bams             // tuple(meta, bam, bai)  meta 含 sample / assembler
    ch_representatives  // tuple(meta, mag_id, mag_fasta)

    main:
    ch_versions = Channel.empty()

    // 列名消歧: 同一样本出现多个 assembler 的 BAM (--assembler both) 时,
    // 列名加 assembler 后缀; 单组装器时列名就是样本 ID
    ch_bam_ambiguous = ch_bams
        .map { meta, bam, bai -> tuple(meta.sample, meta.assembler) }
        .groupTuple()
        .map { sample, assemblers -> tuple(sample, assemblers.unique().size() > 1) }

    // BAM 清单: <BAM basename> \t <列名> —— 进程内列名映射的键 (见 coverm.nf)
    ch_bam_manifest = ch_bams
        .map { meta, bam, bai -> tuple(meta.sample, meta.assembler, bam, bai) }
        .join(ch_bam_ambiguous)
        .map { sample, assembler, bam, bai, ambiguous ->
            def col = ambiguous ? "${sample}.${assembler}" : sample
            "${bam.getFileName()}\t${col}"
        }
        .collectFile(name: 'bam_manifest.tsv', newLine: true, sort: true)

    // 文件依赖: 单 List 传 path 输入 (-resume 依赖追踪 + 任务目录暂存)
    ch_bam_files = ch_bams
        .flatMap { meta, bam, bai -> [bam, bai] }
        .toSortedList()

    ch_mag_files = ch_representatives
        .map { meta, mag_id, mag_fasta -> mag_fasta }
        .toSortedList()

    // 解析脚本作为显式输入暂存 (而非依赖 PATH), 使本 process 可从任意入口复用
    ch_script = file("${projectDir}/bin/coverm_abundance.py", checkIfExists: true)

    // 各组件均为单元素通道, 按位置配对成一次 COVERM 调用; 任一为空
    // (skip 链式传导) 则配对不成立, 整个阶段不调度 (正常收尾)
    COVERM(ch_bam_manifest, ch_bam_files, ch_mag_files, ch_script)
    ch_versions = ch_versions.mix(COVERM.out.versions)

    emit:
    abundance      = COVERM.out.abundance
    multiqc_files  = COVERM.out.abundance  // Phase 15: mag_abundance.tsv 进 MultiQC
    versions       = ch_versions
}
