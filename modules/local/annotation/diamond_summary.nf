// ============================================================================
// DIAMOND_SUMMARY — 解析全部 MAG 的 DIAMOND 命中表为键控汇总表 (Phase 12)
//
// 把本批次全部 <mag_id>.diamond.tsv 经 manifest (meta_id/mag_id/文件映射)
// 与 proteins.faa 序列头校验, 汇总为一张 (meta_id, mag_id, gene) 键控的
// diamond_hits.tsv —— Phase 14 总表 join 的键。gene 列取自 proteins.faa
// 序列头 (Prodigal gene id), 未知 gene id 直接报错, 不静默跳过。
//
// 输入:  bin/diamond_summary.py + manifest + 全部原始命中表 + 全部 faa
// 输出:  diamond_hits.tsv
// ============================================================================

process DIAMOND_SUMMARY {
    label 'process_single'

    conda "conda-forge::python=3.12"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/12_annotation/diamond" }, mode: 'copy'

    input:
    path manifest
    path raw_outs
    path faas
    path script

    output:
    path "diamond_hits.tsv", emit: summary
    path "versions.yml", emit: versions

    script:
    """
    python3 ${script} \\
        --manifest ${manifest} \\
        --output diamond_hits.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """

    stub:
    """
    # stub 也走真实解析脚本 (Phase 10 模式): 上游 stub 的空命中表 + 空 faa
    # 恰好是"空输入"场景, 借此真实验证 manifest 装配、basename 暂存与
    # 输出 schema。键一致性 (gene 取自 faa 头、未知 id 报错) 由 bin/
    # 单测覆盖。
    python3 ${script} \\
        --manifest ${manifest} \\
        --output diamond_hits.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
