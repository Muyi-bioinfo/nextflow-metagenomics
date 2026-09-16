// ============================================================================
// INTEGRATION — 结果整合子工作流 (Phase 14)
//
// 输入:  上游各 Phase 的汇总表通道 (QC/分类/成员表/注释/丰度/代表 MAG/组装
//        汇总)。全部为"只消费"语义 —— 本阶段不重做任何上游计算。
// 处理:  集合级批处理 ——
//        1. 代表 MAG 清单物化为 manifest (collectFile, sort: true), FASTA
//           经 toSortedList 传 path 输入 (ArrayBag 规避, 详见文档);
//        2. 可选表通道以 0 字节哨兵 (assets/empty.tsv) 兜底 —— ifEmpty 只
//           接受具体值, 传通道对象会泄漏 DataflowStream (已知问题), 因此
//           由本子工作流显式提供哨兵文件, 脚本按"0 字节 = 表缺失"处理;
//        3. INTEGRATE_METADATA 产出 mag_metadata.tsv (+ 成员表引用拷贝),
//           INTEGRATE_FUNCTIONAL 产出 mag_functional_annotation.tsv,
//           INTEGRATE_COPIES 拷贝 assembly_summary.tsv。
// 输出:  metadata / functional / membership_copy / assembly_summary /
//        multiqc_files (两张核心表, Phase 15 MultiQC 消费)
//
// 上游 skip 组合导致通道为空时: manifest/FASTA 列表不发射 → 元数据/功能
// 整合不调度 (正常收尾); 拷贝进程独立按各自输入调度。与 Phase 9-13 的
// 空通道行为一致 (上游 skip 已有告警兜底, 不重复)。
// ============================================================================

include { INTEGRATE_METADATA   } from '../../modules/local/integration/metadata.nf'
include { INTEGRATE_FUNCTIONAL } from '../../modules/local/integration/functional.nf'
include { INTEGRATE_COPIES     } from '../../modules/local/integration/copies.nf'

workflow INTEGRATION {

    take:
    ch_representatives
    ch_membership
    ch_qc_table
    ch_taxonomy_table
    ch_abundance
    ch_diamond_table
    ch_eggnog_table
    ch_rgi_table
    ch_assembly_summary

    main:
    ch_versions = Channel.empty()

    // 0 字节哨兵: 可选表的缺省兜底 (脚本按"0 字节 = 表缺失"处理)。
    // ifEmpty 只接受具体值 —— 传通道对象会泄漏 DataflowStream (已知问题),
    // 故显式提供 assets/empty.tsv。
    ch_sentinel = file("${moduleDir}/../../assets/empty.tsv", checkIfExists: true)
    ch_qc_in      = ch_qc_table.ifEmpty(ch_sentinel)
    ch_tax_in     = ch_taxonomy_table.ifEmpty(ch_sentinel)
    ch_abund_in   = ch_abundance.ifEmpty(ch_sentinel)
    ch_diamond_in = ch_diamond_table.ifEmpty(ch_sentinel)
    ch_eggnog_in  = ch_eggnog_table.ifEmpty(ch_sentinel)
    ch_rgi_in     = ch_rgi_table.ifEmpty(ch_sentinel)

    // pathway 映射表: 可选外部依赖 (--pathway_db), 未提供时用哨兵兜底
    // (Pathway 列留空, 不虚构 —— 见 docs/database.md)
    ch_pathway = params.pathway_db
        ? Channel.value(file(params.pathway_db, checkIfExists: true))
        : Channel.value(ch_sentinel)

    // 代表 MAG 清单: meta_id/mag_id/FASTA 绝对路径 (collectFile sort: true
    // → 内容与 -resume 哈希稳定; 绝对路径插值同 Phase 12 manifest 模式)
    ch_rep_manifest = ch_representatives
        .map { meta, mag_id, mag_fasta -> "${meta.id}\t${mag_id}\t${mag_fasta}" }
        .collectFile(name: 'integrate_manifest.tsv', newLine: true, sort: true)

    // FASTA 列表: toSortedList 返回普通 List 而非 ArrayBag (已知问题),
    // 与单元素 manifest 按索引配对
    ch_rep_fastas = ch_representatives
        .map { meta, mag_id, mag_fasta -> mag_fasta }
        .toSortedList()

    ch_metadata_script = file("${moduleDir}/../../bin/integrate_metadata.py", checkIfExists: true)
    INTEGRATE_METADATA(ch_rep_manifest, ch_rep_fastas, ch_membership,
                       ch_qc_in, ch_tax_in, ch_abund_in, ch_metadata_script)
    ch_versions = ch_versions.mix(INTEGRATE_METADATA.out.versions)

    ch_functional_script = file("${moduleDir}/../../bin/integrate_functional.py", checkIfExists: true)
    INTEGRATE_FUNCTIONAL(ch_rep_manifest, ch_diamond_in, ch_eggnog_in,
                         ch_rgi_in, ch_pathway, ch_functional_script)
    ch_versions = ch_versions.mix(INTEGRATE_FUNCTIONAL.out.versions)

    INTEGRATE_COPIES(ch_assembly_summary)

    emit:
    metadata         = INTEGRATE_METADATA.out.metadata
    membership_copy  = INTEGRATE_METADATA.out.membership_copy
    functional       = INTEGRATE_FUNCTIONAL.out.functional
    assembly_summary = INTEGRATE_COPIES.out.assembly_summary
    multiqc_files    = INTEGRATE_METADATA.out.metadata
                           .mix(INTEGRATE_FUNCTIONAL.out.functional)
    versions         = ch_versions
}
