// ============================================================================
// READ_BASED_MERGE — read-based 跨样本整合子工作流 (Phase 20)
//
// 输入:  READ_BASED 的逐样本 emit ——
//        ch_bracken_abundance = tuple(meta, level, path)  (Bracken 各层级丰度表)
//        ch_pathabundance     = tuple(meta, path)         (HUMAnN 通路丰度表)
// 处理:  只消费, 不重做任何上游计算。集合级批处理 ——
//        1. manifest collectFile 物化 (sort: true, 内容与 -resume 哈希稳定),
//           文件列表经 toSortedList 传 path 输入 (ArrayBag 规避, 见 STATUS
//           已知问题; 脚本经 manifest 的 resolve() 双路径打开, Phase 12/14 模式);
//        2. MERGE_BRACKEN 产出各层级 merged_<level>.tsv + beta_diversity.tsv;
//           MERGE_PATHABUNDANCE 产出 merged_pathabundance.tsv。
// 输出:  bracken_merged / beta_diversity / pathabundance_merged / versions
//
// skip 交互: skip_kraken2 / skip_bracken / skip_humann 各自跳过时, 对应逐样本
// emit 为空通道 → collectFile 不发射 → 对应 merge process 不调度 (对应矩阵
// 不产出)。告警已由 read_based.nf 的 skip 守卫发出, 本子工作流不重复。
// ============================================================================

include { MERGE_BRACKEN        } from '../../modules/local/read_based/merge_bracken.nf'
include { MERGE_PATHABUNDANCE  } from '../../modules/local/read_based/merge_pathabundance.nf'

workflow READ_BASED_MERGE {

    take:
    ch_bracken_abundance   // tuple(meta, level, path)
    ch_pathabundance       // tuple(meta, path)

    main:
    ch_versions = Channel.empty()

    // ---------------------------------------------------------------------
    // Bracken 合并: manifest (level \t sample \t 绝对路径) + 文件列表
    // ---------------------------------------------------------------------
    ch_bracken_manifest = ch_bracken_abundance
        .map { meta, level, path -> "${level}\t${meta.id}\t${path}" }
        .collectFile(name: 'bracken_merge_manifest.tsv', newLine: true, sort: true)

    ch_bracken_files = ch_bracken_abundance
        .map { meta, level, path -> path }
        .toSortedList()

    ch_merge_script = file("${moduleDir}/../../bin/merge_read_based.py", checkIfExists: true)
    MERGE_BRACKEN(ch_bracken_manifest, ch_bracken_files, ch_merge_script)
    ch_versions = ch_versions.mix(MERGE_BRACKEN.out.versions)

    // ---------------------------------------------------------------------
    // HUMAnN pathabundance 合并: manifest (sample \t 绝对路径) + 文件列表
    // ---------------------------------------------------------------------
    ch_pathabundance_manifest = ch_pathabundance
        .map { meta, path -> "${meta.id}\t${path}" }
        .collectFile(name: 'pathabundance_merge_manifest.tsv', newLine: true, sort: true)

    ch_pathabundance_files = ch_pathabundance
        .map { meta, path -> path }
        .toSortedList()

    MERGE_PATHABUNDANCE(ch_pathabundance_manifest, ch_pathabundance_files, ch_merge_script)
    ch_versions = ch_versions.mix(MERGE_PATHABUNDANCE.out.versions)

    emit:
    bracken_merged        = MERGE_BRACKEN.out.merged
    beta_diversity        = MERGE_BRACKEN.out.beta_diversity
    pathabundance_merged  = MERGE_PATHABUNDANCE.out.merged
    versions              = ch_versions
}
