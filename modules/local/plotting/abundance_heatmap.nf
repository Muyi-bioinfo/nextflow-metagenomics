// ============================================================================
// PLOT_ABUNDANCE_HEATMAP — MAG 丰度热图 (Phase 21, 图 ①)
//
// 消费 Phase 13 的 mag_abundance.tsv (行 = mag_id, 列 = 样本, 相对丰度 0-1),
// 产出 13_abundance/figures/mag_abundance_heatmap.png。本机真实可验证
// (CoverM 无数据库依赖)。
//
// 只消费现成表, 不重做任何计算。空表 (0 字节 / 仅表头) 时脚本告警跳过,
// PNG 输出为 optional (不产出则 13_abundance/figures/ 不建)。
// ============================================================================

process PLOT_ABUNDANCE_HEATMAP {
    label 'process_single'

    conda "conda-forge::python=3.12 conda-forge::matplotlib=3.10.9 conda-forge::numpy=1.26"
    // 容器: python:3.12 biocontainer 不含 matplotlib/numpy (同 coverm 无 python
    // 的镜像缺口), 容器模式需 mulled 多工具镜像或镜像内自带 matplotlib ——
    // 本机无引擎未实机验证, 待有引擎环境处理 (Phase 17 容器口径)。
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/13_abundance/figures" },
        mode: 'copy', pattern: "*.png"

    input:
    path matrix
    path script

    output:
    path "mag_abundance_heatmap.png", emit: figure, optional: true
    path "versions.yml",              emit: versions

    script:
    """
    python3 ${script} mag-abundance-heatmap \\
        --matrix ${matrix} \\
        --output mag_abundance_heatmap.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """

    stub:
    """
    # stub 走真实脚本 (Phase 20 模式): COVERM stub 的 mag_abundance.tsv 为
    # 0 字节 (按"空表"处理) → 告警跳过, 验证接线与 output schema; 真实 PNG
    # 渲染由 bin/plot_results.py 单测 (合成表) 覆盖, 不伪造真实数值。
    python3 ${script} mag-abundance-heatmap \\
        --matrix ${matrix} \\
        --output mag_abundance_heatmap.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
