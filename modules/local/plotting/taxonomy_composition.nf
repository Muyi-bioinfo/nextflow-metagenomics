// ============================================================================
// PLOT_TAXONOMY_COMPOSITION — MAG 门/纲组成 (Phase 21, 图 ③)
//
// 消费 Phase 10 的 mag_taxonomy.tsv (phylum / class 列), 产出
// 10_mag_taxonomy/figures/mag_taxonomy_composition.png (top-N 门/纲 MAG 计数)。
//
// 只消费现成表。空表/全空分类告警跳过。database-dependent: 真实 GTDB-Tk
// 数据本地无库, 代码由 stub/合成表验证"能画"。
// ============================================================================

process PLOT_TAXONOMY_COMPOSITION {
    label 'process_single'

    conda "conda-forge::python=3.12 conda-forge::matplotlib=3.10.9 conda-forge::numpy=1.26"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/10_mag_taxonomy/figures" },
        mode: 'copy', pattern: "*.png"

    input:
    path taxonomy
    path script

    output:
    path "mag_taxonomy_composition.png", emit: figure, optional: true
    path "versions.yml",                 emit: versions

    script:
    """
    python3 ${script} mag-taxonomy-composition \\
        --taxonomy ${taxonomy} \\
        --output mag_taxonomy_composition.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """

    stub:
    """
    python3 ${script} mag-taxonomy-composition \\
        --taxonomy ${taxonomy} \\
        --output mag_taxonomy_composition.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
