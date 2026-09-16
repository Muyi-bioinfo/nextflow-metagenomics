// ============================================================================
// INTEGRATE_METADATA — 整合 MAG 元数据汇总表 (Phase 14)
//
// 集合级单次调用: 代表 MAG 清单 (manifest) + 成员表 + 各 Phase 汇总表一次
// join 为 mag_metadata.tsv (行 = 代表 MAG, 与 dRep 代表集一致)。可选表
// (QC/分类/丰度) 的缺失由 INTEGRATION 子工作流以 0 字节哨兵文件
// (assets/empty.tsv) 兜底 —— ifEmpty 只接受具体值 (通道对象会泄漏为
// DataflowStream, 详见文档), 脚本按"0 字节 = 表缺失"处理,
// 对应列留空, 与 taxonomy_summary.py 的 --bac/--ar 缺失语义一致。
// Genome_size / GC 由脚本用 Python 标准库现场计算 (不引入新工具)。
//
// 输入:  manifest      = collectFile 物化的 "meta_id \t mag_id \t fasta 绝对
//                        路径" 清单 (sort: true, -resume 哈希稳定)
//        fastas        = 全部代表 MAG FASTA 的单 List (toSortedList, 声明为
//                        path 输入以保证 -resume 依赖追踪与任务目录暂存;
//                        脚本经 manifest 的 resolve() 双路径打开)
//        membership    = Phase 9 成员表 (mag_id \t sample \t rep_mag_id)
//        qc/taxonomy/abundance = 可选 Phase 8/10/13 表 (缺失时为哨兵文件)
//
// 输出:  mag_metadata.tsv
//        mag_membership.tsv  引用拷贝: 成员表原样发布到 14_integrated/
//                            (内容与 09_dereplication 一致, 14_integrated
//                            因此自包含)
// ============================================================================

process INTEGRATE_METADATA {
    tag "mag_metadata"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/14_integrated" }, mode: 'copy'

    input:
    path manifest
    path fastas
    path membership
    // 可选表以 stageAs 固定各自暂存名: 缺失时多个输入兜底为同一个哨兵
    // 文件 (assets/empty.tsv), 若按原名暂存会触发 Nextflow 的输入同名冲突
    path qc_table,       stageAs: 'qc_table.tsv'
    path taxonomy_table, stageAs: 'taxonomy_table.tsv'
    path abundance,      stageAs: 'abundance.tsv'
    path script

    output:
    path "mag_metadata.tsv",  emit: metadata
    path "mag_membership.tsv", emit: membership_copy
    path "versions.yml",       emit: versions

    script:
    """
    python3 ${script} \\
        --manifest ${manifest} \\
        --membership ${membership} \\
        --qc ${qc_table} \\
        --taxonomy ${taxonomy_table} \\
        --abundance ${abundance} \\
        --output mag_metadata.tsv

    # 引用拷贝: 上游 DREP 的成员表产出文件名即 mag_membership.tsv, 暂存后与
    # 本 process 的输出同名 (cp 到自身会报 "same file") —— 文件已在位时直接
    # 作为输出发布 (与 CHECKM2 stub 的同名修复同模式); 名字不同才拷贝。
    if [ -e mag_membership.tsv ]; then
        :
    else
        cp ${membership} mag_membership.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """

    stub:
    """
    # stub 也走真实解析脚本 (Phase 10 模式): 上游 stub 的 qc/分类/成员表均为
    # 真实结构, 丰度 stub 是 0 字节文件 (按"表缺失"处理) —— 借此在 stub-run
    # 真实验证 manifest 装配、resolve() 双路径与输出 schema。Genome_size/GC
    # 的数值正确性由 bin/ 脚本单测覆盖 (stub FASTA 非真实组装)。
    python3 ${script} \\
        --manifest ${manifest} \\
        --membership ${membership} \\
        --qc ${qc_table} \\
        --taxonomy ${taxonomy_table} \\
        --abundance ${abundance} \\
        --output mag_metadata.tsv

    if [ -e mag_membership.tsv ]; then
        :
    else
        cp ${membership} mag_membership.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
