// ============================================================================
// BRACKEN — 基于 Kraken2 报告的丰度估计 (abundance estimation)
//
// ─── 与 Kraken2 的职责分界 ────────────────────────────────────────────────
//   Kraken2  = classification  为每条 read 指派分类单元
//   Bracken  = abundance       重新分配 read 计数, 估计各分类单元的真实丰度
//
// 二者不是"同一件事的两步精修": Kraken2 把大量 read 停留在高层节点 (无法区分
// 同属内的种时, read 被指派到属甚至科)。Bracken 用数据库预先算出的 k-mer 分布,
// 将这些高层 read 按贝叶斯方式重新分配到指定层级, 给出该层级的丰度估计。
//
// 因此 Bracken 只读 Kraken2 的 **report**, 不读逐 read 的 classification output。
//
// ─── 层级 scatter ─────────────────────────────────────────────────────────
// 一次 Bracken 运行只估计一个层级。本模块按 (样本 × 层级) 展开为独立 task
// (见 subworkflows/local/read_based.nf 中的 combine), 因此:
//   - 层级由 params.bracken_levels 配置, 增删层级无需改动模块
//   - 各层级独立缓存, -resume 粒度更细
//
// level 作为独立的 tuple 元素而非并入 meta —— meta 是样本身份的规范表示, 与其他
// 每样本通道 join 时必须保持一致; 把 level 塞进 meta 会让这些 join 失配。
//
// 输入:  tuple(meta, kraken2_report, level)
//        path(db)   Bracken 数据库 (= Kraken2 数据库目录, 内含 databaseXXXmers.kmer_distrib)
// 输出:  tuple(meta, level, abundance_tsv)   该层级的丰度估计表
//        tuple(meta, level, bracken_report)  Kraken 风格报告, 计数已由 Bracken 重估
// ============================================================================

process BRACKEN {
    tag "${meta.id}|${level}"
    label 'process_low'

    // 版本号说明: conda 包为 3.1, 但 `bracken -v` 自报 "Bracken v3.0.1" ——
    // 上游 wrapper 脚本里的 VERSION 字符串未随包版本更新。
    // conda/container 用包版本 (唯一能解析的坐标), versions.yml 记录工具自报值。
    conda "bioconda::bracken=3.1"
    container 'quay.io/biocontainers/bracken:3.1--h9948957_0'

    publishDir { "${params.outdir}/${params.batch_id}/03_taxonomy/bracken" }, mode: 'copy',
        pattern: "*.bracken.*.{tsv,report.txt}"

    input:
    tuple val(meta), path(report), val(level)
    path db

    output:
    tuple val(meta), val(level), path("${meta.id}.bracken.${level}.tsv"),        emit: abundance
    tuple val(meta), val(level), path("${meta.id}.bracken.${level}.report.txt"), emit: report
    path "versions.yml",                                                        emit: versions

    script:
    def read_len = params.bracken_read_length
    """
    # ------------------------------------------------------------------
    # 预检 kmer_distrib —— 必须自己做。
    # bracken 的 wrapper 在数据库检查失败时执行裸 `exit`, 退出码继承上一条
    # echo (即 0), 于是 Nextflow 会把这次失败当成功, 只是产不出文件。
    # 这里显式检查并以非零退出, 让失败真正暴露出来。
    # ------------------------------------------------------------------
    KMER_DISTRIB="${db}/database${read_len}mers.kmer_distrib"

    if [ ! -e "\$KMER_DISTRIB" ]; then
        echo "ERROR: 未找到 Bracken k-mer 分布文件: \$KMER_DISTRIB" >&2
        echo "       Bracken 需要在 Kraken2 数据库上预先构建 (bracken-build -d <db> -l ${read_len})。" >&2
        echo "       文件名中的 ${read_len} 来自 --bracken_read_length, 必须与构建时的 -l 一致;" >&2
        echo "       现有的分布文件为:" >&2
        ls -1 ${db}/*.kmer_distrib 2>/dev/null | sed 's/^/         /' >&2 || echo "         (无)" >&2
        exit 1
    fi

    bracken \\
        -d ${db} \\
        -i ${report} \\
        -o ${meta.id}.bracken.${level}.tsv \\
        -w ${meta.id}.bracken.${level}.report.txt \\
        -r ${read_len} \\
        -l ${level} \\
        -t ${params.bracken_threshold}

    # 同理: bracken 对输入报告的错误也可能以退出码 0 收场, 因此校验产物存在
    if [ ! -s "${meta.id}.bracken.${level}.tsv" ]; then
        echo "ERROR: Bracken 未能生成 ${level} 层级的丰度表 (输出为空或缺失)。" >&2
        echo "       常见原因: Kraken2 报告中该层级没有任何 read。" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: \$( bracken -v 2>&1 | grep -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' | head -n 1 )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.bracken.${level}.tsv ${meta.id}.bracken.${level}.report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: 3.0.1
    END_VERSIONS
    """
}
