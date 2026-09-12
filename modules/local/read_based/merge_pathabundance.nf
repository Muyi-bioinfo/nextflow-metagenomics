// ============================================================================
// MERGE_PATHABUNDANCE — HUMAnN pathabundance → 样本×pathway 宽表矩阵 (Phase 20)
//
// 集合级单次调用: 全部样本的逐样本 pathabundance 表一次合并为一张宽表
// (行 = pathway 字符串, 列 = 样本, 数值 = HUMAnN Abundance (RPK), 缺失
// pathway/样本补 0)。只消费上游 HUMANN 的逐样本 emit, 不重做任何 HUMAnN
// 计算; 本 Phase 只产出数据矩阵, 不画图。
//
// 输入:  manifest = collectFile 物化的 "sample \t 绝对路径" 清单
//                   (sort: true, -resume 哈希稳定)
//        files     = 全部逐样本 pathabundance 表的单 List (toSortedList, 声明为
//                    path 输入以保证 -resume 依赖追踪与任务目录暂存; 脚本经
//                    manifest 的 resolve() 双路径打开, Phase 12/14 模式)
//        script    = bin/merge_read_based.py (显式 path 输入, 不依赖 PATH)
//
// 输出:  merged_pathabundance.tsv
//
// 数据库缺失/跳过时: 上游 HUMANN emit 空通道 → collectFile 不发射 → 本
// process 不调度 (read_based.nf 已告警, 不重复)。
// ============================================================================

process MERGE_PATHABUNDANCE {
    tag "merged_pathabundance"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/04_function/combined" }, mode: 'copy',
        pattern: "*.tsv"

    input:
    path manifest
    path files
    path script

    output:
    path "merged_pathabundance.tsv", emit: merged
    path "versions.yml",             emit: versions

    script:
    """
    python3 ${script} pathabundance \\
        --manifest ${manifest} \\
        --output merged_pathabundance.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """

    stub:
    """
    # stub 也走真实解析脚本 (Phase 12/14 模式): 上游 HUMANN stub 的逐样本表为
    # 0 字节文件 (按"空表"处理), 借此真实验证 manifest 装配、resolve() 双路径
    # 与输出 schema (仅表头宽表)。合并数值正确性由 bin/merge_read_based.py
    # 单测覆盖 (合成表), stub 不伪造真实数据。
    python3 ${script} pathabundance \\
        --manifest ${manifest} \\
        --output merged_pathabundance.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
