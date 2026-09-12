// ============================================================================
// PLOTTING — 结果可视化子工作流 (Phase 21)
//
// 只消费各 Phase 已产出的 TSV (含 Phase 20 combined 矩阵), 产出 7 张 PNG 到
// 对应 Phase 的 figures/ 子目录 (有图才建)。不重做任何上游计算, 不改任何
// 已验证 process —— 仅在下游新增绘图消费层。
//
// 输入 (全部为"只消费"语义):
//   ch_bin_summary         bin_summary.tsv (07)      → 漏斗 raw bins (锚)
//   ch_qc_table            mag_qc.tsv (08)           → 散点 + 漏斗
//   ch_taxonomy_table      mag_taxonomy.tsv (10)     → 组成 + 漏斗
//   ch_membership          mag_membership.tsv (09)   → 漏斗 after dRep
//   ch_abundance           mag_abundance.tsv (13)    → 丰度热图
//   ch_bracken_merged      merged_<level>.tsv (20)   → 跨样本 taxa
//   ch_beta_diversity      beta_diversity.tsv (20, optional) → PCoA
//   ch_pathabundance_merged merged_pathabundance.tsv (20) → pathway 热图
//
// 调度语义: 每个 plot process 以其输入通道是否为空独立调度 —— 对应 Phase 被
// skip 时 (通道空) 该图自然不调度、figures/ 不建。漏斗以 bin_summary 为锚
// (非空才调度), 其余 3 张可选表以 assets/empty.tsv 哨兵兜底 (ifEmpty 只接受
// 具体值, 见 STATUS 已知问题), 脚本按"0 字节 = 层级缺失"跳过该级画剩余漏斗。
// ============================================================================

include { PLOT_ABUNDANCE_HEATMAP    } from '../../modules/local/plotting/abundance_heatmap.nf'
include { PLOT_QC_SCATTER           } from '../../modules/local/plotting/qc_scatter.nf'
include { PLOT_TAXONOMY_COMPOSITION } from '../../modules/local/plotting/taxonomy_composition.nf'
include { PLOT_TAXONOMIC_COMPOSITION } from '../../modules/local/plotting/taxonomic_composition.nf'
include { PLOT_BETA_PCOA            } from '../../modules/local/plotting/beta_pcoa.nf'
include { PLOT_PATHWAY_HEATMAP      } from '../../modules/local/plotting/pathway_heatmap.nf'
include { PLOT_WORKFLOW_SUMMARY     } from '../../modules/local/plotting/workflow_summary.nf'

workflow PLOTTING {

    take:
    ch_bin_summary
    ch_qc_table
    ch_taxonomy_table
    ch_membership
    ch_abundance
    ch_bracken_merged
    ch_beta_diversity
    ch_pathabundance_merged

    main:
    ch_versions = Channel.empty()

    ch_script = file("${projectDir}/bin/plot_results.py", checkIfExists: true)
    ch_sentinel = file("${projectDir}/assets/empty.tsv", checkIfExists: true)

    // ---------------------------------------------------------------------
    // ① MAG 丰度热图 (13_abundance/figures/) —— 本机真实可验证
    // ---------------------------------------------------------------------
    PLOT_ABUNDANCE_HEATMAP(ch_abundance, ch_script)
    ch_versions = ch_versions.mix(PLOT_ABUNDANCE_HEATMAP.out.versions)

    // ---------------------------------------------------------------------
    // ② 完整度 vs 污染度散点 (08_mag_qc/figures/)
    // ---------------------------------------------------------------------
    PLOT_QC_SCATTER(ch_qc_table, ch_script)
    ch_versions = ch_versions.mix(PLOT_QC_SCATTER.out.versions)

    // ---------------------------------------------------------------------
    // ③ MAG 门/纲组成 (10_mag_taxonomy/figures/)
    // ---------------------------------------------------------------------
    PLOT_TAXONOMY_COMPOSITION(ch_taxonomy_table, ch_script)
    ch_versions = ch_versions.mix(PLOT_TAXONOMY_COMPOSITION.out.versions)

    // ---------------------------------------------------------------------
    // ④ 跨样本 top taxa (03_taxonomy/figures/) —— 各层级 merged_*.tsv (glob
    //    emit 已是文件列表), 脚本选 S 前缀层级。直接传通道: 空通道 (skip_
    //    bracken/read_based) 不调度; 单/多文件由 process path 输入 + 脚本
    //    nargs='+' 原样接住 (不用 toSortedList —— 空通道会发空列表误调度)。
    // ---------------------------------------------------------------------
    PLOT_TAXONOMIC_COMPOSITION(ch_bracken_merged, ch_script)
    ch_versions = ch_versions.mix(PLOT_TAXONOMIC_COMPOSITION.out.versions)

    // ---------------------------------------------------------------------
    // ⑤ β 多样性 PCoA (03_taxonomy/figures/) —— optional, 缺失不调度
    // ---------------------------------------------------------------------
    PLOT_BETA_PCOA(ch_beta_diversity, ch_script)
    ch_versions = ch_versions.mix(PLOT_BETA_PCOA.out.versions)

    // ---------------------------------------------------------------------
    // ⑥ 通路丰度热图 (04_function/figures/) —— top-N pathway
    // ---------------------------------------------------------------------
    PLOT_PATHWAY_HEATMAP(ch_pathabundance_merged, ch_script)
    ch_versions = ch_versions.mix(PLOT_PATHWAY_HEATMAP.out.versions)

    // ---------------------------------------------------------------------
    // ⑦ MAG 工作流漏斗 (99_multiqc/figures/) —— bin_summary 为锚, 其余 3 张
    //    可选表以哨兵兜底 (stageAs 固定名规避同名冲突)
    // ---------------------------------------------------------------------
    PLOT_WORKFLOW_SUMMARY(
        ch_bin_summary,
        ch_qc_table.ifEmpty(ch_sentinel),
        ch_membership.ifEmpty(ch_sentinel),
        ch_taxonomy_table.ifEmpty(ch_sentinel),
        ch_script
    )
    ch_versions = ch_versions.mix(PLOT_WORKFLOW_SUMMARY.out.versions)

    emit:
    versions = ch_versions
}
