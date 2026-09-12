// ============================================================================
// PLOT_QC_SCATTER — 完整度 vs 污染度散点 + 阈值线 (Phase 21, 图 ②)
//
// 消费 Phase 8 的 mag_qc.tsv (Completeness / Contamination), 产出
// 08_mag_qc/figures/completeness_vs_contamination.png, 阈值线取
// params.mag_min_completeness / params.mag_max_contamination。
//
// 只消费现成表。空表告警跳过; 单点正常绘制 (不崩)。database-dependent:
// 真实 CheckM2 数据本机无库, 代码由 stub/合成表验证"能画"。
// ============================================================================

process PLOT_QC_SCATTER {
    label 'process_single'

    conda "conda-forge::python=3.12 conda-forge::matplotlib=3.10.9 conda-forge::numpy=1.26"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/08_mag_qc/figures" },
        mode: 'copy', pattern: "*.png"

    input:
    path qc
    path script

    output:
    path "completeness_vs_contamination.png", emit: figure, optional: true
    path "versions.yml",                      emit: versions

    script:
    """
    python3 ${script} qc-scatter \\
        --qc ${qc} \\
        --output completeness_vs_contamination.png \\
        --min-completeness ${params.mag_min_completeness} \\
        --max-contamination ${params.mag_max_contamination}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """

    stub:
    """
    python3 ${script} qc-scatter \\
        --qc ${qc} \\
        --output completeness_vs_contamination.png \\
        --min-completeness ${params.mag_min_completeness} \\
        --max-contamination ${params.mag_max_contamination}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
