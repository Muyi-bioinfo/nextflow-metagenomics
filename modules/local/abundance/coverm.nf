// ============================================================================
// COVERM — MAG 跨样本丰度矩阵 (Phase 13)
//
// genome 模式, 集合级**单次调用**: 一次吃全部样本的 Phase 6 排序 BAM + 全部
// 代表 MAG FASTA, 直接产出矩阵 (行 = MAG, 列 = 样本)。BAM 直接复用 Phase 6
// 产物 —— MAG FASTA 的 contig header 原样保留 (bin/split_bins.py), 与 BAM 的
// contig 名一致, 因此无需重新比对、无需新映射步骤。无外部数据库依赖
// (CoverM 自带比对器)。
//
// 全部逻辑 (BAM 交集过滤 → CoverM 调用 → 键控规范化) 封装在
// bin/coverm_abundance.py 中, 本模块只做薄封装:
//   - BAM 过滤: 与任何 MAG 都无 contig 交集的 BAM (如 0-bin 样本) 会触发
//     CoverM 硬报错, 而非产出全零列 —— 过滤掉, 其样本列由脚本补 0
//   - 规范化: 列 = 样本 ID (both 模式 <sample>.<assembler> 消歧, 映射表来自
//     bam manifest), 行 = mag_id, 剔除 unmapped 伪行, "(%)" 列换算 0-1
// 详见脚本头部注释。
//
// 注意: 单组装器模式下各样本的组装相互独立, contig 名只在同一样本内一致,
// 因此矩阵呈块对角 (MAG 只在自己样本的 BAM 里有非零丰度) —— 这是 single
// 模式的语义而非缺陷, coassembly 接入后矩阵自然变稠密。
//
// ─── 输入 ─────────────────────────────────────────────────────────────────
//   bam_manifest = collectFile 物化的 "<BAM basename> \t <列名>" 清单
//                  (sort: true, -resume 哈希稳定)
//   bam_files    = 全部 BAM + BAI 的单 List (toSortedList, 声明为 path 输入
//                  以保证 -resume 依赖追踪与任务目录暂存)
//   mag_files    = 全部代表 MAG FASTA 的单 List (同上)
//   script       = bin/coverm_abundance.py (checkIfExists, 任何入口可复用)
//
// 注意 coverm_args 不要覆盖 --bam-files / --genome-fasta-files / --methods /
// --output-file / --threads, 它们由脚本固定传入。
// ============================================================================

process COVERM {
    tag "abundance: ${params.coverm_method}"
    label 'process_medium'

    conda "bioconda::coverm=0.8.0 conda-forge::python=3.12"
    // 容器注意: coverm 的 biocontainer 是纯 Rust 产物、不含 python3 —— 脚本
    // 需要 python 解释器, 容器模式请在镜像内提供 python 或换用 mulled 多工具
    // 镜像 (Phase 17 容器验证时处理)
    container 'quay.io/biocontainers/coverm:0.8.0--h750ce8b_0'

    publishDir { "${params.outdir}/${params.batch_id}/13_abundance" }, mode: 'copy'

    input:
    path bam_manifest   // <BAM basename> \t <列名>
    path bam_files      // List: 全部 BAM + BAI
    path mag_files      // List: 全部代表 MAG FASTA (<mag_id>.fa)
    path script

    output:
    path "mag_abundance.tsv", emit: abundance
    path "versions.yml",      emit: versions

    script:
    """
    set -euo pipefail

    python3 ${script} \\
        --bam-manifest ${bam_manifest} \\
        --genomes *.fa \\
        --method ${params.coverm_method} \\
        --threads ${task.cpus} \\
        --coverm-args "${params.coverm_args}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: 0.8.0
    END_VERSIONS
    """

    stub:
    """
    touch mag_abundance.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: 0.8.0
    END_VERSIONS
    """
}
