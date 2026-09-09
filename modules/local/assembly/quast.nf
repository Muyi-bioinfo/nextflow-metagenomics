// ============================================================================
// QUAST — 组装质量评估 (assembly QC)
//
// 职责: 只评估**组装本身**的连续性与规模指标。
//
// ─── 与 CheckM2 的区别 (务必不要混淆) ─────────────────────────────────────
//   QUAST   = assembly QC —— 对象是 contigs 集合, 回答"这次组装拼得好不好":
//             contig 数、总长、N50/N90、最长 contig、GC …
//   CheckM2 = MAG QC     —— 对象是单个分箱基因组, 回答"这个 MAG 是不是一个
//             完整且不混杂的基因组": 完整度 / 污染度 (Phase 8, 与本模块无关)
// N50 高不代表 MAG 质量好; 完整度高也不代表组装连续性好。两者不可互相替代。
//
// 输入:  tuple(meta, contigs)   contigs 可为 .fa 或 .fa.gz (QUAST 原生支持 gz)
// 输出:  results                QUAST 完整输出目录 (report.html / report.tsv /
//                               basic_stats/ / icarus 可视化 …)
//        tsv                    report.tsv 的唯一命名副本, 供跨组装汇总使用
//
// ─── 无参考基因组模式 ─────────────────────────────────────────────────────
// 宏基因组组装没有对应的参考基因组, 因此这里不传 -r。QUAST 退化为纯统计模式,
// 只报告与参考无关的指标 (misassembly / genome fraction 等需要参考的指标不会
// 出现)。metaQUAST 的 --max-ref-number 在线检索参考在生产环境不可靠, 不使用。
//
// ─── N90 与 --report-all-metrics ──────────────────────────────────────────
// QUAST 默认只报告 N50/N75/L50/L75。N90/L90/auN 需要 --report-all-metrics
// (QUAST ≥ 5.2)。该选项在旧版本不存在, 因此这里在运行时探测 --help 后再决定
// 是否附加 —— 不假设容器/conda 解析出的具体版本。
//
// ─── 输出目录命名 ─────────────────────────────────────────────────────────
// 目录名带 ${meta.id}.${meta.assembler} 前缀。这是必需的: MultiQC 靠固定文件名
// report.tsv 发现 QUAST 结果, 若多个组装的目录同名, 暂存时会互相覆盖。
// assembler 进入名字后, --assembler both 的两套结果也能并存。
// ============================================================================

process QUAST {
    tag "${meta.id}.${meta.assembler}"
    label 'process_low'

    conda "bioconda::quast=5.3.0"
    container 'quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2'

    publishDir { "${params.outdir}/${params.batch_id}/05_assembly/quast" },
        mode: 'copy', pattern: "*.quast"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}.${meta.assembler}.quast"),            emit: results
    tuple val(meta), path("${meta.id}.${meta.assembler}.quast.report.tsv"), emit: tsv
    path "versions.yml",                                                    emit: versions

    script:
    def prefix = "${meta.id}.${meta.assembler}"
    """
    # N90/L90/auN 需要 --report-all-metrics (QUAST >= 5.2)。旧版本没有该选项,
    # 直接传会导致参数错误 —— 因此先探测再决定, 缺失时静默降级为默认指标集。
    EXTRA_METRICS=""
    if quast.py --help 2>&1 | grep -q -- "--report-all-metrics"; then
        EXTRA_METRICS="--report-all-metrics"
    else
        echo "WARNING: 当前 QUAST 版本不支持 --report-all-metrics, N90/L90/auN 将缺失。" >&2
    fi

    quast.py \\
        ${contigs} \\
        -o ${prefix}.quast \\
        -l "${prefix}" \\
        --threads ${task.cpus} \\
        --min-contig ${params.quast_min_contig} \\
        \$EXTRA_METRICS \\
        ${params.quast_args}

    # 唯一命名的副本: 跨组装汇总时所有 report.tsv 会被暂存到同一目录,
    # 沿用 QUAST 的固定文件名会互相覆盖。
    cp ${prefix}.quast/report.tsv ${prefix}.quast.report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$( quast.py --version 2>&1 | grep -o 'QUAST v[0-9.]*' | head -n 1 | sed 's/^QUAST v//' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}.${meta.assembler}.quast
    printf 'Assembly\\t${meta.id}.${meta.assembler}\\n# contigs\\t0\\nN50\\t0\\nN90\\t0\\n' \\
        > ${meta.id}.${meta.assembler}.quast/report.tsv
    cp ${meta.id}.${meta.assembler}.quast/report.tsv ${meta.id}.${meta.assembler}.quast.report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: 5.3.0
    END_VERSIONS
    """
}
