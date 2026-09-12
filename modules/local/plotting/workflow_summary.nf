// ============================================================================
// PLOT_WORKFLOW_SUMMARY — MAG 工作流漏斗 (Phase 21, 图 ⑦)
//
// 消费 4 张表数行数, 产出 99_multiqc/figures/workflow_summary.png:
//   raw bins      = bin_summary.tsv 数据行数
//   passing QC    = mag_qc.tsv 中 Completeness >= min 且 Contamination <= max
//   after dRep    = mag_membership.tsv 中 mag_id == representative_mag_id (代表数)
//   GTDB-Tk assigned = mag_taxonomy.tsv 中 domain 非空的行数
// 任一输入层级缺失 (0 字节哨兵 / 0 行) 跳过该级画剩余漏斗; 全缺则整体跳过。
//
// 可选表 (qc/membership/taxonomy) 由 PLOTTING 子工作流以 assets/empty.tsv
// 哨兵兜底, stageAs 固定暂存名规避同名冲突 (同 INTEGRATE_METADATA 模式);
// bin_summary 为锚 (非空才调度本 process —— 无 MAG 流水线则无漏斗)。
// database-dependent: 真实 CheckM2/GTDB-Tk 数据本机无库。
// ============================================================================

process PLOT_WORKFLOW_SUMMARY {
    label 'process_single'

    conda "conda-forge::python=3.12 conda-forge::matplotlib=3.10.9 conda-forge::numpy=1.26"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/99_multiqc/figures" },
        mode: 'copy', pattern: "*.png"

    input:
    path bin_summary
    path qc_table,       stageAs: 'qc_table.tsv'
    path membership,     stageAs: 'membership.tsv'
    path taxonomy,       stageAs: 'taxonomy.tsv'
    path script

    output:
    path "workflow_summary.png", emit: figure, optional: true
    path "versions.yml",          emit: versions

    script:
    """
    python3 ${script} workflow-summary \\
        --bin-summary ${bin_summary} \\
        --qc ${qc_table} \\
        --membership ${membership} \\
        --taxonomy ${taxonomy} \\
        --output workflow_summary.png \\
        --min-completeness ${params.mag_min_completeness} \\
        --max-contamination ${params.mag_max_contamination}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """

    stub:
    """
    python3 ${script} workflow-summary \\
        --bin-summary ${bin_summary} \\
        --qc ${qc_table} \\
        --membership ${membership} \\
        --taxonomy ${taxonomy} \\
        --output workflow_summary.png \\
        --min-completeness ${params.mag_min_completeness} \\
        --max-contamination ${params.mag_max_contamination}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
