// ============================================================================
// DEREPLICATION — MAG 去冗余子工作流 (Phase 9)
//
// (工作流名为 DEREPLICATION, 与 process DREP 区分 —— 同名会在模块内冲突;
// mag.nf 以 `DEREPLICATION as DREP` 引入, 接线为 DREP.out.*)
//
// 输入:  tuple(meta, mag_id, mag_fasta)  (全部 qualified MAG)
//        path(qc_table)                  (Phase 8 的 mag_qc.tsv)
//
// 处理:  集合级批处理 —— 全部 qualified MAG 聚合成一个基因组目录, dRep 一次
//        dereplicate 聚类为代表集合; 解析聚类输出生成成员表。
// 输出:  representatives = tuple(meta, mag_id, mag_fasta)   形状同输入,
//                                                            仅代表 MAG
//        membership      = path(mag_membership.tsv)          mag_id \t sample
//                                                            \t representative
//        catalog         = path(dereplicated_genomes/)       代表 MAG FASTA 目录
//        clusters        = path(data_tables/)                dRep 原始聚类数据表
//
// 聚合规则 (ArrayBag 规避, 详见文档):
//   - manifest/genomeInfo 类元数据: map{...}.collectFile(name:..., sort: true)
//     物化为文件;
//   - 文件列表: map{...}.toSortedList() 传 path 输入。
// 空 qualified 通道时 collectFile/toSortedList 都不发射, DREP/DREP_SUMMARY
// 不调度, representatives/membership 均为空 —— 与 Phase 10 空通道行为一致
// (上游 skip_binning 已有告警兜底)。
//
// --skip_dereplication: representatives 直接传递 qualified_mags, membership
// 生成恒等映射 (每个 MAG 代表自身), 输出等价于全量分类。
// ============================================================================

include { DREP         } from '../../modules/local/dereplication/drep.nf'
include { DREP_SUMMARY } from '../../modules/local/dereplication/derep_summary.nf'

workflow DEREPLICATION {

    take:
    ch_qualified_mags
    ch_qc_table

    main:
    ch_versions = Channel.empty()

    if (params.skip_dereplication) {
        log.warn "MAG 去冗余已跳过 (--skip_dereplication), 全部 qualified MAG 直接作为代表 MAG, 成员表为恒等映射 (等价全量分类)。"

        ch_representatives = ch_qualified_mags
        ch_membership = ch_qualified_mags
            .map { meta, mag_id, mag_fasta -> "${mag_id}\t${meta.id}\t${mag_id}" }
            .collectFile(name: 'mag_membership.tsv', newLine: true, sort: true)
        ch_catalog  = Channel.empty()
        ch_clusters = Channel.empty()
    }
    else {
        // dRep 选代表按 completeness/contamination 打分, 质量输入来自 Phase 8
        // 的 QC 表; --skip_mag_qc 时没有 QC 表可构建 genomeInfo (dRep 会退化为
        // 任意选代表), 明确报错而非静默降级
        if (params.skip_mag_qc) {
            error "ERROR: MAG dereplication (dRep) requires CheckM2 QC results to build genomeInfo. Either provide --checkm2_db (do not skip MAG QC) or pass --skip_dereplication."
        }

        // manifest: mag_id → sample (成员表的 sample 列; 同时是解析脚本的
        // 输入)。sort: true 保证文件内容与 -resume 哈希稳定。
        ch_manifest = ch_qualified_mags
            .map { meta, mag_id, mag_fasta -> "${mag_id}\t${meta.id}" }
            .collectFile(name: 'genome_manifest.tsv', newLine: true, sort: true)

        // genomeInfo.csv (genome,completeness,contamination): 从 QC 表映射,
        // 且只保留 qualified MAG 的行 (QC 表含全部输入 MAG, 含被阈值过滤的)。
        // dRep 的 genome 列不带 .fa 扩展名, genome 即 mag_id。
        ch_genome_info = ch_qc_table
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.mag_id, "${row.mag_id},${row.Completeness},${row.Contamination}") }
            .join(ch_qualified_mags.map { meta, mag_id, mag_fasta -> tuple(mag_id) })
            .map { mag_id, genome_info_row -> genome_info_row }
            .collectFile(name: 'genomeInfo.csv', newLine: true, sort: true)

        // MAG FASTA 路径列表 (toSortedList 返回普通 List 而非 ArrayBag)
        ch_fastas = ch_qualified_mags
            .map { meta, mag_id, mag_fasta -> mag_fasta }
            .toSortedList()

        DREP(ch_genome_info, ch_fastas)
        ch_catalog  = DREP.out.catalog
        ch_clusters = DREP.out.clusters
        ch_versions = ch_versions.mix(DREP.out.versions)

        ch_parser = file("${projectDir}/bin/derep_summary.py", checkIfExists: true)
        DREP_SUMMARY(ch_manifest, DREP.out.clusters, DREP.out.catalog, ch_parser)
        ch_membership = DREP_SUMMARY.out.membership
        ch_versions = ch_versions.mix(DREP_SUMMARY.out.versions)

        // 代表 MAG 恢复: 成员表中 mag_id == representative_mag_id 的行即代表,
        // 与原始 qualified 通道按 mag_id join 找回 meta 与 FASTA —— 形状保持
        // tuple(meta, mag_id, mag_fasta), Phase 10 直接消费
        ch_rep_ids = DREP_SUMMARY.out.membership
            .splitCsv(header: false, sep: '\t')
            .filter { row -> row[0] == row[2] }
            .map { row -> tuple(row[0]) }

        ch_representatives = ch_qualified_mags
            .map { meta, mag_id, mag_fasta -> tuple(mag_id, meta, mag_fasta) }
            .join(ch_rep_ids)
            .map { mag_id, meta, mag_fasta -> tuple(meta, mag_id, mag_fasta) }
    }

    emit:
    representatives = ch_representatives
    membership      = ch_membership
    catalog         = ch_catalog
    clusters        = ch_clusters
    versions        = ch_versions
}
