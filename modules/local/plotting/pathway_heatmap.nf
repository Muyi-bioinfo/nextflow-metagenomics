// ============================================================================
// PLOT_PATHWAY_HEATMAP — 通路丰度热图 (Phase 21, 图 ⑥)
//
// 消费 Phase 20 的 merged_pathabundance.tsv (样本×pathway 宽表), 取 top-N
// pathway (N = params.plot_pathway_top, 默认 50) 按丰度排序, 产出
// 04_function/figures/pathway_abundance_heatmap.png (log10 色阶)。
//
// 只消费现成表。空表/无 pathway 告警跳过。database-dependent: 真实 HUMAnN
// 数据本地无库, 代码由 stub/合成表验证"能画"。
// ============================================================================

process PLOT_PATHWAY_HEATMAP {
    label 'process_single'

    conda "conda-forge::python=3.12 conda-forge::matplotlib=3.10.9 conda-forge::numpy=1.26"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/04_function/figures" },
        mode: 'copy', pattern: "*.png"

    input:
    path matrix
    path script

    output:
    path "pathway_abundance_heatmap.png", emit: figure, optional: true
    path "versions.yml",                  emit: versions

    script:
    """
    python3 ${script} pathway-heatmap \\
        --matrix ${matrix} \\
        --output pathway_abundance_heatmap.png \\
        --top ${params.plot_pathway_top}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """

    stub:
    """
    python3 ${script} pathway-heatmap \\
        --matrix ${matrix} \\
        --output pathway_abundance_heatmap.png \\
        --top ${params.plot_pathway_top}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
