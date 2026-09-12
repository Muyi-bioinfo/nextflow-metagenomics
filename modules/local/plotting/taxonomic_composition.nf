// ============================================================================
// PLOT_TAXONOMIC_COMPOSITION — 跨样本 top taxa (Phase 21, 图 ④)
//
// 消费 Phase 20 的 merged_<level>.tsv (样本×taxa 宽表, 各层级一张), 脚本选
// S 前缀层级 (否则首个), 产出 03_taxonomy/figures/taxonomic_composition.png
// (top-N taxa 按跨样本总丰度排序)。
//
// 只消费现成表。空表/无 taxa 告警跳过。database-dependent: 真实 Bracken 数据
// 本机无库, 代码由 stub/合成表验证"能画"。
// ============================================================================

process PLOT_TAXONOMIC_COMPOSITION {
    label 'process_single'

    conda "conda-forge::python=3.12 conda-forge::matplotlib=3.10.9 conda-forge::numpy=1.26"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/03_taxonomy/figures" },
        mode: 'copy', pattern: "*.png"

    input:
    path matrices    // merged_<level>.tsv 文件列表 (各层级各一张)
    path script

    output:
    path "taxonomic_composition.png", emit: figure, optional: true
    path "versions.yml",               emit: versions

    script:
    """
    python3 ${script} taxonomic-composition \\
        --matrices ${matrices} \\
        --output taxonomic_composition.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """

    stub:
    """
    python3 ${script} taxonomic-composition \\
        --matrices ${matrices} \\
        --output taxonomic_composition.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
