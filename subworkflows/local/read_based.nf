// ============================================================================
// READ_BASED — 基于 reads 的分析子工作流 (无需组装)
//
// ─── 数据流: 两条并行分支 ─────────────────────────────────────────────────
//
//   clean_reads ──┬─→ KRAKEN2 ──→ report ──→ BRACKEN (× 各分类层级)
//                 │      分类              丰度估计
//                 │
//                 └─→ HUMANN
//                        功能谱
//
// 两条分支同时从 `clean_reads` 取数据, 彼此没有依赖 —— Nextflow 会并行调度。
//
// 明确不是这样:
//   KRAKEN2 → BRACKEN → HUMANN        (串联)
// HUMAnN 不消费 Kraken2/Bracken 的任何产物: 它自带 MetaPhlAn prescreen 决定
// 检索哪些物种的基因组。把它挂在 Bracken 之后只会白白串行化, 且 Kraken2 报告
// 与 HUMAnN 的 --taxonomic-profile 格式不兼容, 根本传不进去。
//
// ─── 分支内部的真实依赖 ───────────────────────────────────────────────────
// BRACKEN 在 KRAKEN2 之后, 这是数据依赖而非设计缺陷: Bracken 的输入就是
// Kraken2 的 report。职责划分见各模块头注释:
//   KRAKEN2 = classification (逐 read 指派分类单元)
//   BRACKEN = abundance      (把停留在高层节点的 read 重新分配到目标层级)
//
// ─── 层级 scatter ─────────────────────────────────────────────────────────
// 一次 Bracken 运行只估计一个层级。这里用 combine 把 (样本) × (层级) 展开成
// 独立 task —— params.bracken_levels 默认 "S,G" 即种 + 属两级丰度。
//
// ─── 数据库缺失时的行为 ───────────────────────────────────────────────────
// 每个工具各自检查自己的数据库: 缺失则 log.warn 并跳过该分支, 其余分支照常运行。
// 绝不伪造产物 —— 跳过的分支 emit 空通道, 日志中有明确原因。
// (若用户显式给了路径但路径不存在, 则由 checkIfExists 硬失败 —— 那是配置错误。)
//
// ─── 输出目录 ─────────────────────────────────────────────────────────────
//   03_taxonomy/kraken2/   Kraken2 报告 (+ 可选的逐 read 分类结果)
//   03_taxonomy/bracken/   各层级丰度估计表与 Bracken 报告
//   04_function/humann/    基因家族 / 通路丰度 / 通路覆盖度
// ============================================================================

include { KRAKEN2 } from '../../modules/local/read_based/kraken2.nf'
include { BRACKEN } from '../../modules/local/read_based/bracken.nf'
include { HUMANN  } from '../../modules/local/read_based/humann.nf'

workflow READ_BASED {

    take:
    ch_reads    // channel: [ val(meta), path(reads) ]  reads=[R1,R2] (PE) 或 [R1] (SE)

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // 跳过的分支 emit 空通道 —— emit 块中的每个通道都必须存在, 与是否执行无关
    ch_kraken2_report   = Channel.empty()
    ch_kraken2_output   = Channel.empty()
    ch_bracken_abundance = Channel.empty()
    ch_bracken_report   = Channel.empty()
    ch_genefamilies     = Channel.empty()
    ch_pathabundance    = Channel.empty()
    ch_pathcoverage     = Channel.empty()

    // =====================================================================
    // 分支 A: 物种分类与丰度 (Kraken2 → Bracken)
    // =====================================================================
    def run_kraken2 = params.kraken2_db && !params.skip_kraken2

    if (run_kraken2) {
        // 数据库目录作为 value channel: 每个样本复用同一次暂存。
        // 用 path 输入而非字符串, 让 Nextflow 以符号链接暂存 (数据库可达数百 GB)。
        ch_kraken2_db = Channel.value(file(params.kraken2_db, checkIfExists: true))

        KRAKEN2(ch_reads, ch_kraken2_db)

        ch_kraken2_report = KRAKEN2.out.report
        ch_kraken2_output = KRAKEN2.out.classified_output
        ch_versions       = ch_versions.mix(KRAKEN2.out.versions.first())

        // MultiQC 的 kraken 模块靠报告内容 (而非文件名) 识别。
        // 注意只送 Kraken2 的报告: Bracken 的报告是同一格式但计数已重估,
        // 一并送入会被当成额外样本, 造成重复计数。
        ch_multiqc_files = ch_multiqc_files.mix(
            KRAKEN2.out.report.map { meta, report -> report }
        )

        // -----------------------------------------------------------------
        // Bracken —— 输入是 Kraken2 的 report, 因此必然在 Kraken2 之后
        // -----------------------------------------------------------------
        // kmer_distrib 文件通常就构建在 Kraken2 数据库目录内, 故默认复用该路径
        def bracken_db = params.bracken_db ?: params.kraken2_db

        if (!params.skip_bracken) {
            // 层级列表解析与校验: 拼错 (如写成 "species") 会在 task 深处才报错,
            // 这里提前拦下。Bracken 接受主层级首字母 + 可选子层级数字 (如 S1)。
            def levels = params.bracken_levels.toString()
                .tokenize(',')
                .collect { lvl -> lvl.trim() }
                .findAll { lvl -> lvl }

            if (!levels) {
                error "ERROR: --bracken_levels 未指定任何层级 (当前值: '${params.bracken_levels}')。"
            }

            def bad_levels = levels.findAll { lvl -> !(lvl ==~ /^[RKDPCOFGS][0-9]*$/) }
            if (bad_levels) {
                error "ERROR: --bracken_levels 含无效层级 ${bad_levels}。" +
                      " 应使用 Kraken 层级代码 R/K/D/P/C/O/F/G/S, 可带子层级数字 (如 S1)。" +
                      " 例: --bracken_levels 'S,G' (种 + 属)"
            }

            // scatter: (样本 × 层级) → [ meta, report, level ]
            // level 保持为独立元素而非并入 meta —— meta 是样本身份的规范表示,
            // 与其他每样本通道 join 时必须一致。
            ch_bracken_input = ch_kraken2_report.combine(Channel.fromList(levels))

            BRACKEN(ch_bracken_input, Channel.value(file(bracken_db, checkIfExists: true)))

            ch_bracken_abundance = BRACKEN.out.abundance
            ch_bracken_report    = BRACKEN.out.report
            ch_versions          = ch_versions.mix(BRACKEN.out.versions.first())
        }
        else {
            log.warn "Bracken 已跳过 (--skip_bracken)。仅有 Kraken2 的 classification 结果, 没有丰度估计。"
        }
    }
    else {
        if (params.skip_kraken2) {
            log.warn "Kraken2 已跳过 (--skip_kraken2), 不产出物种分类结果。"
        }
        else {
            log.warn "未提供 --kraken2_db, Kraken2 已跳过, 不产出物种分类结果。"
        }

        // Bracken 的输入是 Kraken2 报告, 因此必须一并跳过 —— 显式说明级联原因,
        // 避免用户以为 --bracken_db 已配好就会有丰度结果。
        if (!params.skip_bracken) {
            log.warn "Bracken 随之跳过 —— 它的输入是 Kraken2 报告, 无 Kraken2 结果即无法进行丰度估计。"
        }
    }

    // =====================================================================
    // 分支 B: 功能谱 (HUMAnN) —— 与分支 A 并行, 不依赖其任何产物
    // =====================================================================
    // HUMAnN 需要三个数据库:
    //   nucleotide (ChocoPhlAn) / protein (UniRef) / MetaPhlAn (prescreen 用)
    // CLAUDE.md 预留的单个 humann_db 不足以表达 —— 这里保留它作为"父目录"便捷
    // 写法, 按 `humann_databases --download` 的标准布局推导子目录, 且仅在子目录
    // 确实存在时才采用 (不猜测路径); 三个专用参数可逐个覆盖。
    def resolve_humann_subdir = { String subdir ->
        if (!params.humann_db) return null
        def candidate = file("${params.humann_db}/${subdir}")
        return candidate.exists() ? candidate.toString() : null
    }

    def humann_nt_db  = params.humann_nucleotide_db ?: resolve_humann_subdir.call('chocophlan')
    def humann_pr_db  = params.humann_protein_db    ?: resolve_humann_subdir.call('uniref')
    def humann_mpa_db = params.metaphlan_db

    if (!params.skip_humann && humann_nt_db && humann_pr_db && humann_mpa_db) {
        HUMANN(
            ch_reads,
            Channel.value(file(humann_nt_db,  checkIfExists: true)),
            Channel.value(file(humann_pr_db,  checkIfExists: true)),
            Channel.value(file(humann_mpa_db, checkIfExists: true))
        )

        ch_genefamilies  = HUMANN.out.genefamilies
        ch_pathabundance = HUMANN.out.pathabundance
        ch_pathcoverage  = HUMANN.out.pathcoverage
        ch_versions      = ch_versions.mix(HUMANN.out.versions.first())

        // MultiQC 没有 HUMAnN 模块, 因此不向 ch_multiqc_files 追加任何内容。
    }
    else if (params.skip_humann) {
        log.warn "HUMAnN 已跳过 (--skip_humann), 不产出功能谱结果。"
    }
    else {
        // 逐个点明缺哪个, 而不是笼统说"数据库不全"。
        // 单行拼接而非多行: Nextflow 的控制台输出会把换行折叠成一行。
        def missing = []
        if (!humann_nt_db)  missing << "nucleotide/ChocoPhlAn (--humann_nucleotide_db 或 --humann_db/chocophlan)"
        if (!humann_pr_db)  missing << "protein/UniRef (--humann_protein_db 或 --humann_db/uniref)"
        if (!humann_mpa_db) missing << "MetaPhlAn (--metaphlan_db)"

        log.warn "HUMAnN 已跳过 —— 缺少数据库: " + missing.join(' | ') +
                 " 。不产出功能谱结果 (基因家族 / 通路丰度 / 通路覆盖度)。"
    }

    emit:
    // ---- 分支 A: 03_taxonomy ----
    kraken2_report    = ch_kraken2_report      // [ val(meta), path(report) ]
    kraken2_output    = ch_kraken2_output      // [ val(meta), path(output.txt.gz) ]
    bracken_abundance = ch_bracken_abundance   // [ val(meta), val(level), path(tsv) ]
    bracken_report    = ch_bracken_report      // [ val(meta), val(level), path(report) ]

    // ---- 分支 B: 04_function ----
    genefamilies      = ch_genefamilies        // [ val(meta), path(tsv) ]
    pathabundance     = ch_pathabundance       // [ val(meta), path(tsv) ]
    pathcoverage      = ch_pathcoverage        // [ val(meta), path(tsv) ]

    // ---- 汇总 ----
    multiqc_files     = ch_multiqc_files       // path(*)  仅 Kraken2 报告
    versions          = ch_versions            // path(versions.yml)
}
