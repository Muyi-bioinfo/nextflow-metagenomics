// ============================================================================
// MAG_QC — MAG 质控子工作流 (Phase 8)
//
// 输入:  tuple(meta, mag_id, mag_fasta)
// 输出:  qualified_mags (同形状, 按阈值过滤) + qc_table (mag_qc.tsv)
// ============================================================================

include { CHECKM2         } from '../../modules/local/mag_qc/checkm2.nf'
include { MAG_QC_SUMMARY  } from '../../modules/local/mag_qc/mag_qc_summary.nf'

workflow MAG_QC {

    take:
    ch_mags

    main:
    ch_versions = Channel.empty()

    if (params.skip_mag_qc) {
        log.warn "MAG 质控已跳过 (--skip_mag_qc), qualified_mags 直接传递输入 MAG。"
        ch_qualified = ch_mags
        ch_qc_table = Channel.empty()
    }
    else {
        if (!params.checkm2_db) {
            error "ERROR: --checkm2_db is required when MAG QC is enabled. Provide the CheckM2 DIAMOND database path."
        }

        CHECKM2(ch_mags)
        ch_versions = ch_versions.mix(CHECKM2.out.versions)

        // 汇总为单一 QC 表。collect() 保留空通道语义，由 summary 输出表头。
        MAG_QC_SUMMARY(CHECKM2.out.results.map { meta, mag_id, mag_fasta, qc -> qc }.collect())
        ch_qc_table = MAG_QC_SUMMARY.out.tsv
        ch_versions = ch_versions.mix(MAG_QC_SUMMARY.out.versions)

        // QC 结果与原始 MAG 按 meta.id + mag_id 配对，保留原始 MAG FASTA。
        ch_qc_by_key = CHECKM2.out.results
            .map { meta, mag_id, mag_fasta, qc -> tuple("${meta.id}::${mag_id}", meta, mag_id, mag_fasta, qc) }
        ch_qualified = ch_qc_by_key
            .filter { key, meta, mag_id, mag_fasta, qc ->
                def row = qc.text.readLines().drop(1)[0].split('\\t', -1)
                def completeness = row[3] as BigDecimal
                def contamination = row[4] as BigDecimal
                completeness >= params.mag_min_completeness && contamination <= params.mag_max_contamination
            }
            .map { key, meta, mag_id, mag_fasta, qc -> tuple(meta, mag_id, mag_fasta) }
    }

    emit:
    qualified_mags = ch_qualified
    qc_table       = ch_qc_table
    multiqc_files  = ch_qc_table          // Phase 15: mag_qc.tsv 进 MultiQC
    versions       = ch_versions
}
