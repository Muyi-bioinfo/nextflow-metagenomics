// ============================================================================
// FASTP_SUMMARY — 汇总各样本 fastp 指标为单一 TSV
//
// 职责: 将逐样本的 fastp JSON 展平为一张表, 便于下游消费与人工审阅。
//       实际解析逻辑在 bin/parse_fastp_json.py。
//
// 记录指标: raw/filtered reads 与 bases、保留率、Q20/Q30、GC、duplication、
//           adapter 修剪量、各类过滤失败计数。
//
// 本表是 Phase 14 `sample_qc_summary.tsv` 的组成部分。
// ============================================================================

process FASTP_SUMMARY {
    label 'process_single'

    conda "conda-forge::python=3.12"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/01_qc" }, mode: 'copy'

    input:
    path json_files
    // 解析脚本作为显式输入暂存, 而非依赖 PATH 或 ${projectDir}:
    // 二者都只在「从仓库根调用」时成立, 子工作流被其他入口复用时会失效。
    path parser

    output:
    path "fastp_summary.tsv", emit: tsv

    script:
    """
    python3 ${parser} ${json_files} --output fastp_summary.tsv
    """

    stub:
    """
    touch fastp_summary.tsv
    """
}
