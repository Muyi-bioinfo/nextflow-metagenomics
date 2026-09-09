// ============================================================================
// ASSEMBLY — 从头组装子工作流
//
// ─── 数据流 ───────────────────────────────────────────────────────────────
//
//   clean_reads ──┬─→ MEGAHIT ────→ contigs ──┐
//                 │                            ├─→ QUAST ─→ ASSEMBLY_SUMMARY
//                 └─→ METASPADES ─→ contigs ──┘   组装QC      跨组装汇总表
//
// ─── 关键设计: 两个组装器是「可选替代」, 不是流水线的两级 ─────────────────
//
// 明确不是这样:
//     MEGAHIT → metaSPAdes            (串联)
//
// 而是由 params.assembler 选择:
//     "megahit"     只跑 MEGAHIT      (默认)
//     "metaspades"  只跑 metaSPAdes
//     "both"        两者**并行**各跑一次, 得到两套独立的 contigs
//
// "both" 不是"先 A 后 B", 两条分支同时从 clean_reads 取数据, 互不消费对方产物。
// 它的用途是在同一份 reads 上横向比较两个组装器 —— assembly_summary.tsv 里同一
// 样本会有两行, 可直接对比 N50 / 总长 / contig 数。
//
// ─── 组装单元 (assembly unit) 与 assembly_mode ────────────────────────────
//
// 「组装单元」= 一次组装运行的输入集合。它由 assembly_mode 决定:
//
//   single      每个样本各自组装          S01 → assembly_S01
//                                        S02 → assembly_S02
//   coassembly  多个样本合并后一次组装     S01 + S02 + S03 → coassembly_01
//
// 本 Phase 只实现 single。为了让 coassembly 后续能平滑接入, 组装单元的身份
// 被显式写进 meta 而非隐含在样本 id 里:
//
//   meta.assembler      产出该 contigs 的组装器 (megahit | metaspades)
//   meta.assembly_mode  single | coassembly
//   meta.samples        参与本次组装的样本 id 列表 —— single 时长度为 1,
//                       coassembly 时为该组的全部样本
//
// meta.samples 是关键: Phase 6 需要把 reads 比回 contigs, 必须知道「这套 contigs
// 由哪些样本的 reads 组装而来」。single 模式下它等价于 [meta.id], 但下游代码只
// 依赖 meta.samples, 因此 coassembly 接入时下游无需改动。
//
// coassembly 未实现的部分是**通道分组与 reads 合并**, 不是模块:
//   1. 按 meta.group (或新增的 coassembly_group 字段) groupTuple
//   2. 合并每组的 R1/R2 —— MEGAHIT 接受逗号分隔多文件, metaSPAdes 需先 cat
//   3. 为组装单元生成稳定 id (如 "coassembly_<group>")
// 在此之前显式报错, 而不是静默按 single 处理。
//
// ─── 通道结构 ─────────────────────────────────────────────────────────────
//
//   tuple(meta, contigs)
//
// 而非 tuple(meta, assembly_id, contigs): assembly_id 会与 meta 里的身份信息
// 重复, 且下游 join 时必须两边都记得带上它。组装身份放进 meta 后, 整条流水线
// 的通道形状保持一致 (与 Phase 3/4 的 tuple(meta, ...) 相同), join 直接以 meta
// 为键即可。--assembler both 时同一样本的两个 meta 因 assembler 字段不同而天然
// 可区分, 不会互相覆盖。
//
// ─── 输出目录 ─────────────────────────────────────────────────────────────
//   05_assembly/megahit/<sample>/     contigs + 运行日志
//   05_assembly/metaspades/<sample>/  contigs + scaffolds + 日志 (+ 可选 GFA)
//   05_assembly/quast/<id>.<asm>.quast/   QUAST 完整报告
//   05_assembly/assembly_summary.tsv      跨组装汇总表
// ============================================================================

include { MEGAHIT          } from '../../modules/local/assembly/megahit.nf'
include { METASPADES       } from '../../modules/local/assembly/metaspades.nf'
include { QUAST            } from '../../modules/local/assembly/quast.nf'
include { ASSEMBLY_SUMMARY } from '../../modules/local/assembly/assembly_summary.nf'

// ---------------------------------------------------------------------------
// 由样本 reads 通道构建某个组装器的组装单元通道。
//
// single 模式下是 1:1 映射, 看似多余的一层封装 —— 但它是 coassembly 的唯一接入
// 点: 未来只需在此按组 groupTuple 并合并 reads, 调用方与下游都不必改动。
// ---------------------------------------------------------------------------
def build_assembly_units(ch_reads, String assembler) {
    return ch_reads.map { meta, reads ->
        def assembly_meta = meta + [
            assembler:     assembler,
            assembly_mode: 'single',
            samples:       [ meta.id ]
        ]
        tuple(assembly_meta, reads)
    }
}

workflow ASSEMBLY {

    take:
    ch_reads    // channel: [ val(meta), path(reads) ]  reads=[R1,R2] (PE) 或 [R1] (SE)

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // 跳过的分支 emit 空通道 —— emit 块中的每个通道都必须存在, 与是否执行无关
    ch_contigs          = Channel.empty()
    ch_scaffolds        = Channel.empty()
    ch_quast_results    = Channel.empty()
    ch_assembly_summary = Channel.empty()

    // =====================================================================
    // 参数校验 —— 全部在启动时完成, 不留到 task 深处才失败
    // =====================================================================
    def assembler = params.assembler?.toString()?.toLowerCase()?.trim()
    def valid_assemblers = ['megahit', 'metaspades', 'both']

    if (!valid_assemblers.contains(assembler)) {
        error "ERROR: --assembler 取值无效 (当前: '${params.assembler}')。" +
              " 可选: ${valid_assemblers.join(' | ')}。" +
              " 注意 megahit 与 metaspades 是**可选替代**组装器, 'both' 表示并行各跑一次以作比较。"
    }

    def mode = params.assembly_mode?.toString()?.toLowerCase()?.trim()

    if (mode == 'coassembly') {
        error "ERROR: --assembly_mode coassembly 尚未实现 (当前 Phase 仅实现 single)。" +
              " 通道结构已按可扩展方式设计 (meta.assembly_mode / meta.samples), 但样本分组与 reads 合并逻辑未完成。" +
              " 此处显式报错而非静默按 single 处理 —— 后者会产出与预期完全不同的结果。"
    }

    if (mode != 'single') {
        error "ERROR: --assembly_mode 取值无效 (当前: '${params.assembly_mode}')。可选: single | coassembly(未实现)。"
    }

    // MEGAHIT 的 --presets 本身就是一组预设 k, 与 --k-list 互斥; 同时给会被 MEGAHIT 拒绝
    if (params.megahit_k_list && params.megahit_preset) {
        error "ERROR: --megahit_k_list 与 --megahit_preset 互斥 (preset 本身即一组预设 k 值)。请只指定其中之一。"
    }

    def run_megahit    = assembler in ['megahit', 'both']
    def run_metaspades = assembler in ['metaspades', 'both']

    if (assembler == 'both') {
        log.info "组装器: MEGAHIT 与 metaSPAdes 将**并行**各跑一次 (--assembler both), 产出两套独立 contigs 供比较, 而非串联执行。"
    }

    // =====================================================================
    // 分支 A: MEGAHIT
    // =====================================================================
    if (run_megahit) {
        MEGAHIT(build_assembly_units(ch_reads, 'megahit'))

        ch_contigs  = ch_contigs.mix(MEGAHIT.out.contigs)
        ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())
    }

    // =====================================================================
    // 分支 B: metaSPAdes —— 与分支 A 并行, 不消费其任何产物
    // =====================================================================
    if (run_metaspades) {
        // metaSPAdes (--meta) 只接受双端 library, 这是工具的硬性限制。
        // 单端样本在此分流并明确告知, 而不是让 SPAdes 在 task 里报晦涩的错。
        ch_metaspades_units = build_assembly_units(ch_reads, 'metaspades')
            .branch { meta, reads ->
                paired: !meta.single_end
                single: meta.single_end
            }

        ch_metaspades_units.single.subscribe { meta, reads ->
            log.warn "metaSPAdes 跳过单端样本 ${meta.id} —— metaSPAdes (--meta) 要求双端 library。该样本请改用 --assembler megahit。"
        }

        METASPADES(ch_metaspades_units.paired)

        ch_contigs   = ch_contigs.mix(METASPADES.out.contigs)
        ch_scaffolds = METASPADES.out.scaffolds
        ch_versions  = ch_versions.mix(METASPADES.out.versions.first())
    }

    // =====================================================================
    // 组装 QC: QUAST
    //
    // 每个组装单元一个 QUAST 任务 —— 两个组装器的结果都流经这里, 因此
    // --assembler both 时同一样本会得到两份可直接对比的 QUAST 报告。
    //
    // 再次强调: QUAST 评估的是 contigs 集合的连续性 (N50/N90/...), 不是 MAG
    // 质量。MAG 的完整度/污染度由 Phase 8 的 CheckM2 负责。
    // =====================================================================
    if (!params.skip_quast) {
        QUAST(ch_contigs)

        ch_quast_results = QUAST.out.results
        ch_versions      = ch_versions.mix(QUAST.out.versions.first())

        // MultiQC 的 quast 模块靠固定文件名 report.tsv 发现结果, 因此这里送入
        // 整个 QUAST 输出目录 (目录名含 id + assembler, 暂存时不会互相覆盖)。
        ch_multiqc_files = ch_multiqc_files.mix(
            QUAST.out.results.map { meta, results -> results }
        )

        // 跨组装汇总表 (跨样本 + 跨组装器, 因此 collect)
        // 解析脚本以 moduleDir 相对定位到本仓库, 与调用入口的 projectDir 无关
        ch_quast_parser = file("${moduleDir}/../../bin/parse_quast_report.py", checkIfExists: true)
        ASSEMBLY_SUMMARY(QUAST.out.tsv.map { meta, tsv -> tsv }.collect(), ch_quast_parser)

        ch_assembly_summary = ASSEMBLY_SUMMARY.out.tsv
    }
    else {
        log.warn "QUAST 已跳过 (--skip_quast), 不产出组装质量指标 (contig 数 / 总长 / N50 / N90 / 最长 contig / GC), assembly_summary.tsv 也不会生成。"
    }

    emit:
    // 下游 (Phase 6 比对 → Phase 7 分箱) 唯一的 contigs 来源
    // [ val(meta), path(contigs.fa.gz) ]
    // meta 额外携带: assembler / assembly_mode / samples
    contigs          = ch_contigs

    // metaSPAdes 的 scaffolds —— 供人工审阅, 不进入分箱数据流 (N gap 会干扰
    // 四核苷酸频率与覆盖度计算)
    scaffolds        = ch_scaffolds

    // ---- 组装 QC ----
    quast_results    = ch_quast_results      // [ val(meta), path(quast_dir) ]
    assembly_summary = ch_assembly_summary   // path(assembly_summary.tsv)

    // ---- 汇总 ----
    multiqc_files    = ch_multiqc_files      // path(*)  QUAST 输出目录
    versions         = ch_versions           // path(versions.yml)
}
