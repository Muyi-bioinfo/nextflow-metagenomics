// ============================================================================
// INTEGRATE_FUNCTIONAL — 整合功能注释汇总表 (Phase 14)
//
// 集合级单次调用: 三张 (meta_id, mag_id, gene) 键控注释表按外连接合并为
// 一张 mag_functional_annotation.tsv (行 = 三表键并集, 跨表同键合并为一行,
// 缺失表/缺失字段留空)。可选表的缺失以 0 字节哨兵文件兜底 (与
// INTEGRATE_METADATA 同机制); pathway_db 缺失时 Pathway 列留空 (不伪造)。
// manifest 用作 MAG 身份校验 (未知 mag_id 报错) 与三表全缺时的调度触发
// (输出仅表头, 与 Phase 12 的空输入语义一致)。
//
// 输入:  manifest   = 代表 MAG 清单 (collectFile, 同 INTEGRATE_METADATA)
//        diamond/eggnog/rgi = 可选 Phase 12 注释表 (缺失时为哨兵文件)
//        pathway     = 可选 KO→pathway 映射 (params.pathway_db, 缺失时为
//                      哨兵文件)
//
// 输出:  mag_functional_annotation.tsv
//        (列: MAG_ID / Gene / KO / COG / GO / Pathway / ARG)
// ============================================================================

process INTEGRATE_FUNCTIONAL {
    tag "mag_functional_annotation"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/14_integrated" }, mode: 'copy'

    input:
    path manifest
    // 可选表以 stageAs 固定各自暂存名: 缺失时多个输入兜底为同一个哨兵
    // 文件 (assets/empty.tsv), 若按原名暂存会触发 Nextflow 的输入同名冲突
    path diamond_table, stageAs: 'diamond_table.tsv'
    path eggnog_table,  stageAs: 'eggnog_table.tsv'
    path rgi_table,     stageAs: 'rgi_table.tsv'
    path pathway_db,    stageAs: 'pathway_db.tsv'
    path script

    output:
    path "mag_functional_annotation.tsv", emit: functional
    path "versions.yml", emit: versions

    script:
    """
    python3 ${script} \\
        --manifest ${manifest} \\
        --diamond ${diamond_table} \\
        --eggnog ${eggnog_table} \\
        --rgi ${rgi_table} \\
        --pathway ${pathway_db} \\
        --output mag_functional_annotation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """

    stub:
    """
    # stub 也走真实解析脚本 (Phase 10 模式): 上游 stub 的注释表为真实结构
    # (空命中即"空输入"场景), 借此真实验证三表外连接、键并集与输出 schema。
    # 键一致性 (未知 gene/mag_id 报错) 由 bin/ 脚本单测覆盖。
    python3 ${script} \\
        --manifest ${manifest} \\
        --diamond ${diamond_table} \\
        --eggnog ${eggnog_table} \\
        --rgi ${rgi_table} \\
        --pathway ${pathway_db} \\
        --output mag_functional_annotation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
