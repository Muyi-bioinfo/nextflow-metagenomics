// ============================================================================
// ANNOTATION — Phase 12 功能注释 (DIAMOND blastp ∥ eggNOG-mapper ∥ RGI)
//
// 三支彼此独立的并行分支, 唯一输入都是 GENE_PREDICTION.out.proteins
// (tuple(meta, mag_id, proteins.faa), 每个代表 MAG 一条)。逐 MAG 直接
// 映射, 无集合级聚合。每支产出:
//   原始工具输出 (<mag_id>.<tool>.*, 发布保留) + 一张解析汇总表
//   (meta_id, mag_id, gene 键控, gene 取自 proteins.faa 序列头 ——
//   Phase 14 总表 join 的键, 键不一致即本阶段验收失败)。
//
// 参数守卫 (Phase 4/8/10/11 模式): 数据库为 null 且未加对应 skip 时明确
// 报错, 不静默跳过。skip_annotation 总开关; skip_diamond / skip_eggnog /
// skip_rgi 分开关, skip 时告警并跳过对应分支, 汇总表通道为空。
//
// 上游通道为空 (如 --skip_gene_prediction) 时各分支自动不调度任务
// (collectFile 空通道不发射, 汇总 process 连带不调度), 正常收尾。
// ============================================================================

include { DIAMOND_BLASTP  } from '../../modules/local/annotation/diamond.nf'
include { DIAMOND_SUMMARY } from '../../modules/local/annotation/diamond_summary.nf'
include { EGGNOG_MAPPER  } from '../../modules/local/annotation/eggnog.nf'
include { EGGNOG_SUMMARY } from '../../modules/local/annotation/eggnog_summary.nf'
include { RGI_LOAD; RGI_MAIN } from '../../modules/local/annotation/rgi.nf'
include { RGI_SUMMARY     } from '../../modules/local/annotation/rgi_summary.nf'

workflow ANNOTATION {
    take:
    ch_proteins    // tuple(meta, mag_id, proteins.faa)

    main:
    ch_versions      = Channel.empty()
    ch_diamond_table = Channel.empty()
    ch_eggnog_table  = Channel.empty()
    ch_rgi_table     = Channel.empty()

    if (!params.skip_annotation) {

        // -----------------------------------------------------------------
        // DIAMOND blastp 分支 (NR 蛋白库)
        // -----------------------------------------------------------------
        if (params.skip_diamond) {
            log.warn "DIAMOND 注释已跳过 (--skip_diamond), 不产出 diamond_hits.tsv。"
        }
        else {
            if (!params.diamond_db) {
                error "ERROR: DIAMOND 注释需要 --diamond_db (NR 蛋白库的 .dmnd 文件); 请提供数据库, 或 --skip_diamond 跳过该分支。"
            }
            DIAMOND_BLASTP(ch_proteins)
            DIAMOND_SUMMARY(
                // manifest 写绝对路径: 真实运行中 raw/faa 作为输入暂存进任务
                // 目录, 解析脚本优先按暂存名 (basename) 打开; -stub-run 不
                // 暂存未被 stub 引用的输入, 此时回退到绝对路径 (上游任务
                // work 目录) —— 两种模式解析脚本都能取到文件。
                DIAMOND_BLASTP.out.hits
                    .map { meta, mag_id, out, faa -> "${meta.id}\t${mag_id}\t${out}\t${faa}" }
                    .collectFile(name: 'diamond_manifest.tsv', sort: true, newLine: true),
                // 文件通道必须 toSortedList 聚合成单元素 List: 与单元素
                // manifest 按索引配对, 否则只取到第一个 MAG 的文件 (哪个
                // MAG 取决于任务完成顺序, -resume 哈希不稳定)。排序亦保证
                // 哈希确定性。
                DIAMOND_BLASTP.out.hits.map { meta, mag_id, out, faa -> out }.toSortedList(),
                DIAMOND_BLASTP.out.hits.map { meta, mag_id, out, faa -> faa }.toSortedList(),
                file("${projectDir}/bin/diamond_summary.py", checkIfExists: true)
            )
            ch_diamond_table = DIAMOND_SUMMARY.out.summary
            ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions, DIAMOND_SUMMARY.out.versions)
        }

        // -----------------------------------------------------------------
        // eggNOG-mapper 分支 (COG/GO/EC/KO/Pathway)
        // -----------------------------------------------------------------
        if (params.skip_eggnog) {
            log.warn "eggNOG-mapper 注释已跳过 (--skip_eggnog), 不产出 eggnog_annotations.tsv。"
        }
        else {
            if (!params.eggnog_db) {
                error "ERROR: eggNOG-mapper 注释需要 --eggnog_db (eggNOG 数据目录, 含 eggnog.db 与 eggnog_proteins.dmnd); 请提供数据库, 或 --skip_eggnog 跳过该分支。"
            }
            EGGNOG_MAPPER(ch_proteins)
            EGGNOG_SUMMARY(
                EGGNOG_MAPPER.out.annotations
                    .map { meta, mag_id, ann, faa -> "${meta.id}\t${mag_id}\t${ann}\t${faa}" }
                    .collectFile(name: 'eggnog_manifest.tsv', sort: true, newLine: true),
                // toSortedList 聚合: 见 DIAMOND 分支注释 (单元素 manifest
                // 配对 + 哈希确定性)。
                EGGNOG_MAPPER.out.annotations.map { meta, mag_id, ann, faa -> ann }.toSortedList(),
                EGGNOG_MAPPER.out.annotations.map { meta, mag_id, ann, faa -> faa }.toSortedList(),
                file("${projectDir}/bin/eggnog_summary.py", checkIfExists: true)
            )
            ch_eggnog_table = EGGNOG_SUMMARY.out.summary
            ch_versions = ch_versions.mix(EGGNOG_MAPPER.out.versions, EGGNOG_SUMMARY.out.versions)
        }

        // -----------------------------------------------------------------
        // RGI (CARD) 分支 (ARG/ARO)
        // -----------------------------------------------------------------
        if (params.skip_rgi) {
            log.warn "RGI (CARD) 注释已跳过 (--skip_rgi), 不产出 rgi_annotations.tsv。"
        }
        else {
            if (!params.card_db) {
                error "ERROR: RGI (CARD) 注释需要 --card_db (card.json 路径); 请提供数据库, 或 --skip_rgi 跳过该分支。"
            }
            // RGI_LOAD 以 proteins 首个元素作调度门: 上游为空时不载库。
            // RGI_LOAD.out.db 单元素经 combine 广播给每个 MAG 的 RGI_MAIN
            // (队列通道单元素只能被消费一次, 不能直接与逐 MAG 通道配对)。
            RGI_LOAD(ch_proteins.first())
            RGI_MAIN(ch_proteins.combine(RGI_LOAD.out.db))
            RGI_SUMMARY(
                RGI_MAIN.out.annotations
                    .map { meta, mag_id, rgi_json, faa -> "${meta.id}\t${mag_id}\t${rgi_json}\t${faa}" }
                    .collectFile(name: 'rgi_manifest.tsv', sort: true, newLine: true),
                // toSortedList 聚合: 见 DIAMOND 分支注释 (单元素 manifest
                // 配对 + 哈希确定性)。
                RGI_MAIN.out.annotations.map { meta, mag_id, rgi_json, faa -> rgi_json }.toSortedList(),
                RGI_MAIN.out.annotations.map { meta, mag_id, rgi_json, faa -> faa }.toSortedList(),
                file("${projectDir}/bin/rgi_summary.py", checkIfExists: true)
            )
            ch_rgi_table = RGI_SUMMARY.out.summary
            ch_versions = ch_versions.mix(RGI_LOAD.out.versions, RGI_MAIN.out.versions, RGI_SUMMARY.out.versions)
        }
    }
    else {
        log.warn "整个功能注释阶段已跳过 (--skip_annotation), 不产出 DIAMOND / eggNOG / RGI 注释。"
    }

    emit:
    diamond_table = ch_diamond_table
    eggnog_table  = ch_eggnog_table
    rgi_table     = ch_rgi_table
    versions      = ch_versions
}
