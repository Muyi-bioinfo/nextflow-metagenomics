// ============================================================================
// EGGNOG_SUMMARY — 解析全部 MAG 的 emapper 注释为键控汇总表 (Phase 12)
//
// 把本批次全部 <mag_id>.emapper.annotations 经 manifest (meta_id/mag_id/
// 文件映射) 与 proteins.faa 序列头校验, 汇总为一张 (meta_id, mag_id,
// gene) 键控的 eggnog_annotations.tsv —— Phase 14 总表 join 的键, 携带
// COG_category / Description / Preferred_name / EC / KO / GO /
// KEGG_Pathway (列按列名定位, 版本列序变动不破坏解析)。
//
// 输入:  bin/eggnog_summary.py + manifest + 全部原始注释 + 全部 faa
// 输出:  eggnog_annotations.tsv
// ============================================================================

process EGGNOG_SUMMARY {
    label 'process_single'

    conda "conda-forge::python=3.12"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/12_annotation/eggnog" }, mode: 'copy'

    input:
    path manifest
    path raw_outs
    path faas
    path script

    output:
    path "eggnog_annotations.tsv", emit: summary
    path "versions.yml", emit: versions

    script:
    """
    python3 ${script} \\
        --manifest ${manifest} \\
        --output eggnog_annotations.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """

    stub:
    """
    # stub 也走真实解析脚本 (Phase 10 模式): 上游 stub 的空 annotations +
    # 空 faa 恰好是"空输入"场景, 借此真实验证 manifest 装配、basename
    # 暂存与输出 schema。键一致性由 bin/ 单测覆盖。
    python3 ${script} \\
        --manifest ${manifest} \\
        --output eggnog_annotations.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
