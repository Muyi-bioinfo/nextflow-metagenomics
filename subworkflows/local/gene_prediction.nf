// ============================================================================
// GENE_PREDICTION — 基因预测子工作流 (Phase 11)
//
// 输入:  ch_mags  tuple(meta, mag_id, mag_fasta) —— 来自 DREP.out.representatives,
//                 每个代表 MAG 一条记录。不用 DREP.out.catalog: 它丢 meta,
//                 且 --skip_dereplication 时为空通道 (representatives 此时
//                 恒等传递全部 qualified MAG, 永不空, 天然免兜底)。
// 输出:  proteins = tuple(meta, mag_id, proteins.faa)  ← Phase 12 三支注释的唯一输入
//        genes    = tuple(meta, mag_id, genes.fna)
//        gff      = tuple(meta, mag_id, gff)
//
// Prodigal 无外部数据库依赖, 逐 MAG 运行 (-p meta); 上游通道已是
// tuple(meta, mag_id, mag_fasta) 形状, 直接映射, 无需集合级聚合。
// ============================================================================

include { PRODIGAL } from '../../modules/local/gene_prediction/prodigal.nf'

workflow GENE_PREDICTION {

    take:
    ch_mags

    main:
    ch_versions = Channel.empty()

    if (params.skip_gene_prediction) {
        log.warn "基因预测已跳过 (--skip_gene_prediction), 不产出 genes.fna / proteins.faa / GFF, Phase 12 功能注释也无从进行。"
        ch_proteins = Channel.empty()
        ch_genes    = Channel.empty()
        ch_gff      = Channel.empty()
    }
    else {
        PRODIGAL(ch_mags)
        ch_proteins = PRODIGAL.out.proteins
        ch_genes    = PRODIGAL.out.genes
        ch_gff      = PRODIGAL.out.gff
        ch_versions = ch_versions.mix(PRODIGAL.out.versions)
    }

    emit:
    proteins = ch_proteins
    genes    = ch_genes
    gff      = ch_gff
    versions = ch_versions
}
