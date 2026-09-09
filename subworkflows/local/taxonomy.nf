// ============================================================================
// TAXONOMY — MAG 分类子工作流 (Phase 10)
//
// 输入:  ch_mags      待分类 MAG: tuple(meta, mag_id, mag_fasta) —— 执行顺序
//                    按最佳实践为先去冗余后分类, 故正常路径下是 Phase 9 dRep
//                    产出的**代表 MAG** 集合。
//        ch_membership  dRep 成员表文件 (mag_id \t sample \t representative_
//                    mag_id, 覆盖全部 qualified MAG), **必须由调用方提供**。
//                    Phase 9 接入前, workflows/mag.nf 用恒等映射 (每个 MAG
//                    代表自身) 临时充当, 接入后换成 DREP.out.membership。
// 输出:  mag_taxonomy.tsv (sample, mag_id, rep_mag_id, domain..species) ——
//        覆盖全部 qualified MAG: 代表 MAG 用自身分类, 冗余成员 MAG 的分类
//        由同簇代表 MAG 回填 (rep_mag_id 列标明来源)。
//
// GTDB-Tk classify_wf 是批处理工具: 把输入 MAG 聚合为一个基因组目录, 一次
// 调用完成分类; 聚合时以 MAG ID 命名 FASTA, 使 summary 的 user_genome 列可
// 映射回 (meta.id, mag_id), 并在 manifest 中显式记录该映射。
// ============================================================================

include { GTDBTK           } from '../../modules/local/taxonomy/gtdbtk.nf'
include { TAXONOMY_SUMMARY } from '../../modules/local/taxonomy/taxonomy_summary.nf'

workflow TAXONOMY {

    take:
    ch_mags
    ch_membership

    main:
    ch_versions = Channel.empty()

    if (params.skip_taxonomy) {
        log.warn "MAG 分类已跳过 (--skip_taxonomy), 不产出 mag_taxonomy.tsv。"
        ch_taxonomy_table = Channel.empty()
    }
    else {
        if (!params.gtdbtk_db) {
            error "ERROR: --gtdbtk_db is required when MAG taxonomy is enabled. Provide the GTDB-Tk reference data path (R220+, ~110 GB unpacked)."
        }

        ch_script = file("${projectDir}/bin/taxonomy_summary.py", checkIfExists: true)

        // 聚合: manifest (mag_id ↔ meta.id) 经 collectFile 物化为文件, 排序
        // 保证 -resume 哈希稳定; FASTA 路径经 toSortedList 聚合成单次批处理
        // 输入。空通道时两者都不发射, GTDBTK / TAXONOMY_SUMMARY 均不调度
        // (上游已就无 MAG 告警)。
        // 注意不用 collect() 后传 val(tuple 列表) —— Nextflow 会把聚合结果
        // 扁平化为 ArrayBag, process 内无法按三元组还原 (见 gtdbtk.nf 头注)。
        ch_manifest = ch_mags
            .map { meta, mag_id, mag_fasta -> "${mag_id}\t${meta.id}" }
            .collectFile(name: 'genome_manifest.tsv', newLine: true, sort: true)
        ch_fastas = ch_mags
            .map { meta, mag_id, mag_fasta -> mag_fasta }
            .toSortedList()

        GTDBTK(ch_manifest, ch_fastas)
        ch_versions = ch_versions.mix(GTDBTK.out.versions)

        TAXONOMY_SUMMARY(ch_script, ch_manifest, ch_membership, GTDBTK.out.bac_summary, GTDBTK.out.ar_summary)
        ch_taxonomy_table = TAXONOMY_SUMMARY.out.tsv
        ch_versions = ch_versions.mix(TAXONOMY_SUMMARY.out.versions)
    }

    emit:
    taxonomy_table = ch_taxonomy_table
    multiqc_files  = ch_taxonomy_table    // Phase 15: mag_taxonomy.tsv 进 MultiQC
    versions       = ch_versions
}
