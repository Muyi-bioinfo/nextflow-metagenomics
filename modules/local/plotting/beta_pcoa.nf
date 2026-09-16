// ============================================================================
// PLOT_BETA_PCOA — β 多样性 PCoA (Phase 21, 图 ⑤)
//
// 消费 Phase 20 的 beta_diversity.tsv (样本×样本 Bray-Curtis 距离矩阵),
// numpy 做 Gower 双中心化 + np.linalg.eigh 特征分解, 取前 2 主坐标,
// 产出 03_taxonomy/figures/beta_diversity_pcoa.png (轴标解释方差比例)。
//
// 只消费现成表。距离矩阵缺失 / <2 样本 / 正特征值 <2 时脚本告警跳过。
// beta_diversity 本身是 optional emit (单样本/0 taxa 不发射) → 通道空则本
// process 不调度。database-dependent: 真实 Bracken 数据本地无库。
// ============================================================================

process PLOT_BETA_PCOA {
    label 'process_single'

    conda "conda-forge::python=3.12 conda-forge::matplotlib=3.10.9 conda-forge::numpy=1.26"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/03_taxonomy/figures" },
        mode: 'copy', pattern: "*.png"

    input:
    path distance
    path script

    output:
    path "beta_diversity_pcoa.png", emit: figure, optional: true
    path "versions.yml",             emit: versions

    script:
    """
    python3 ${script} beta-pcoa \\
        --distance ${distance} \\
        --output beta_diversity_pcoa.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """

    stub:
    """
    python3 ${script} beta-pcoa \\
        --distance ${distance} \\
        --output beta_diversity_pcoa.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
