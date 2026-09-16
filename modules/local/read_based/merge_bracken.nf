// ============================================================================
// MERGE_BRACKEN — Bracken 各层级丰度表 → 样本×taxa 宽表矩阵 (Phase 20)
//
// 集合级单次调用: 全部样本 × 全部层级的逐样本 Bracken 丰度表一次合并。
// 每个层级产出一张宽表 (行 = 分类单元名, 列 = 样本, 数值 = fraction_total_reads
// 相对丰度, 缺失样本补 0), 并基于所选层级 (S 前缀优先, 否则首个层级) 的
// 相对丰度向量算 Bray-Curtis 距离矩阵 (样本×样本)。
//
// 只消费上游 BRACKEN 的逐样本 emit, 不重做任何 Bracken 计算; 本 Phase 只产出
// 数据矩阵, 不画图 (画图另行一个 Phase, 消费现成表)。
//
// 输入:  manifest = collectFile 物化的 "level \t sample \t 绝对路径" 清单
//                   (sort: true, -resume 哈希稳定)
//        files     = 全部逐样本 Bracken 丰度表的单 List (toSortedList, 声明为
//                    path 输入以保证 -resume 依赖追踪与任务目录暂存; 脚本经
//                    manifest 的 resolve() 双路径打开, Phase 12/14 模式)
//        script    = bin/merge_read_based.py (显式 path 输入, 不依赖 PATH)
//
// 输出:  merged_<level>.tsv  每个层级一张 (glob 匹配, emit merged)
//        beta_diversity.tsv  Bray-Curtis 距离矩阵 (单样本/0 taxa 时跳过,
//                            optional 输出不产出)
//
// 数据库缺失/跳过时: 上游 BRACKEN emit 空通道 → collectFile 不发射 → 本
// process 不调度 (read_based.nf 已告警, 不重复)。
// ============================================================================

process MERGE_BRACKEN {
    tag "merged_bracken"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/03_taxonomy/combined" }, mode: 'copy',
        pattern: "*.tsv"

    input:
    path manifest
    path files
    path script

    output:
    path "merged_*.tsv",          emit: merged
    path "beta_diversity.tsv",    emit: beta_diversity, optional: true
    path "versions.yml",          emit: versions

    script:
    """
    python3 ${script} bracken \\
        --manifest ${manifest} \\
        --output-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """

    stub:
    """
    # stub 也走真实解析脚本 (Phase 12/14 模式): 上游 BRACKEN stub 的逐样本表为
    # 0 字节文件 (按"空表"处理), 借此真实验证 manifest 装配、resolve() 双路径
    # 与输出 schema (仅表头宽表 + beta diversity 跳过)。合并数值正确性由
    # bin/merge_read_based.py 单测覆盖 (合成表), stub 不虚构真实数据。
    python3 ${script} bracken \\
        --manifest ${manifest} \\
        --output-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
