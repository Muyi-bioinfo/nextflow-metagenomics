// ============================================================================
// TAXONOMY_SUMMARY — 解析 GTDB-Tk summary 为 per-MAG 分类表 (Phase 10)
//
// 把 gtdbtk.bac120.summary.tsv / gtdbtk.ar53.summary.tsv 的 user_genome 列
// 经 genome_manifest.tsv 映射回 (meta.id, mag_id), 拆分 classification 列为
// domain..species 7 级 (缺级为空, -1 split), 并按 dRep 成员表把代表 MAG 的
// 分类回填给同簇冗余成员 MAG (rep_mag_id 列标明来源)。两个 summary 输入恒
// 存在 (GTDBTK process 保证), 缺其一只是空文件, 解析脚本自行跳过。
//
// 输入:  bin/taxonomy_summary.py + manifest + membership + bac/ar summary
// 输出:  mag_taxonomy.tsv
// ============================================================================

process TAXONOMY_SUMMARY {
    label 'process_single'

    conda "conda-forge::python=3.12"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/10_mag_taxonomy" }, mode: 'copy'

    input:
    path script
    path manifest
    path membership
    path bac_summary
    path ar_summary

    output:
    path "mag_taxonomy.tsv", emit: tsv
    path "versions.yml", emit: versions

    script:
    """
    python3 ${script} \\
        --manifest ${manifest} \\
        --output mag_taxonomy.tsv \\
        --membership ${membership} \\
        --bac ${bac_summary} \\
        --ar ${ar_summary}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """

    stub:
    """
    # stub 也走真实解析脚本: 输入 (manifest + 成员表 + 真假 summary) 已就位,
    # 借此在 stub-run 中真实验证 user_genome → (meta.id, mag_id) 的映射与
    # 成员分类回填逻辑。
    python3 ${script} \\
        --manifest ${manifest} \\
        --output mag_taxonomy.tsv \\
        --membership ${membership} \\
        --bac ${bac_summary} \\
        --ar ${ar_summary}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
