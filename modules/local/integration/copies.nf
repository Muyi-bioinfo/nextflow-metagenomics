// ============================================================================
// INTEGRATE_COPIES — 样本级表引用拷贝 (Phase 14)
//
// 把 Phase 5 的 assembly_summary.tsv 原样发布到 14_integrated/。样本级表
// (QUAST/Kraken2/Bracken/HUMAnN) 已在各自阶段产出并发布, 不在 MAG 表 join
// 范围 —— 14_integrated 只做引用拷贝, 不重做。纯文件拷贝, 无外部工具依赖。
//
// --skip_quast / --skip_assembly 时上游通道为空, 本 process 不调度, 自然
// 不产生拷贝 (与"无可拷贝内容"语义一致)。
//
// 输入:  assembly_summary  (Phase 5 QUAST 汇总, 单元素通道或空)
// 输出:  assembly_summary.tsv
// ============================================================================

process INTEGRATE_COPIES {
    tag "assembly_summary"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/14_integrated" }, mode: 'copy'

    input:
    path assembly_summary

    output:
    path "assembly_summary.tsv", emit: assembly_summary

    script:
    """
    # 上游 ASSEMBLY_SUMMARY 的产出文件名即 assembly_summary.tsv, 暂存后与
    # 本 process 的输出同名 (cp 到自身会报 "same file") —— 文件已在位时
    # 直接作为输出发布 (与 CHECKM2 stub 的同名修复同模式); 名字不同才拷贝。
    if [ -e assembly_summary.tsv ]; then
        :
    else
        cp ${assembly_summary} assembly_summary.tsv
    fi
    """

    stub:
    """
    if [ -e assembly_summary.tsv ]; then
        :
    else
        cp ${assembly_summary} assembly_summary.tsv
    fi
    """
}
