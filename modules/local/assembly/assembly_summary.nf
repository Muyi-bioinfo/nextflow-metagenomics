// ============================================================================
// ASSEMBLY_SUMMARY — 汇总各组装的 QUAST 指标为单一 TSV
//
// 职责: 把逐组装的 QUAST report.tsv 展平成一张跨样本 (跨组装器) 的表。
//       实际解析逻辑在 bin/parse_quast_report.py。
//
// 与 FASTP_SUMMARY 同一模式: 单个工具的报告适合人读, 跨样本比较需要一张表。
// --assembler both 时同一样本会有两行 (megahit / metaspades), 便于直接比较
// 两个组装器在同一份 reads 上的表现 —— 这正是 "both" 模式存在的意义。
//
// 记录指标: contig 数、总长、最长 contig、N50/N75/N90、L50/L75/L90、auN、
//           GC%、各长度阈值以上的 contig 数与碱基数、每 100 kbp 的 N 数。
//
// 本表是 Phase 14 `assembly_summary.tsv` 的来源。
// ============================================================================

process ASSEMBLY_SUMMARY {
    label 'process_single'

    conda "conda-forge::python=3.12"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/05_assembly" }, mode: 'copy'

    input:
    path report_files
    // 解析脚本作为显式输入暂存, 而非依赖 PATH 或 ${projectDir}:
    // 二者都只在「从仓库根调用」时成立, 子工作流被其他入口复用时会失效。
    path parser

    output:
    path "assembly_summary.tsv", emit: tsv

    script:
    """
    python3 ${parser} ${report_files} --output assembly_summary.tsv
    """

    stub:
    """
    touch assembly_summary.tsv
    """
}
