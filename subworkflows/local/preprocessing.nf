// ============================================================================
// PREPROCESSING — 预处理子工作流
//
// 数据流:
//   raw_reads ──┬─→ FastQC ────────────────────────────→ (QC 报告)
//               └─→ fastp → clean_reads → HOST_REMOVAL → nonhost_reads
//                              │                              │
//                              └──────────┬───────────────────┘
//                                         ↓
//                                  emit: reads (单一来源)
//
// ─── SE/PE 设计 ───────────────────────────────────────────────────────────
// 通道统一使用 tuple(meta, reads), reads 为列表:
//   PE: [R1, R2]   SE: [R1]
// meta.single_end (Boolean) 由 check_samplesheet.py 写入, 各模块凭此决定命令行。
//
// ─── 关键设计: clean reads 单一来源 ────────────────────────────────────────
// 本子工作流只 emit 一个 `reads` 通道。read-based 分支 (Kraken2/Bracken/HUMAnN)
// 与 assembly 分支 (MEGAHIT/...) 都从这一个通道取数据 —— Nextflow 通道可被多个
// 下游 process 消费, 不会重复计算上游任务。因此两条路线在物理上使用完全相同的
// FASTQ 文件, 而非各自重新生成一套。
//
// ─── 宿主去除的条件执行 ───────────────────────────────────────────────────
// 未提供 --host_index (或显式 --skip_host_removal) 时, HOST_REMOVAL 被跳过,
// `reads` 直接由 fastp 输出承接。此时 emit 的语义仍然成立 (clean reads),
// 但未经宿主去除 —— 日志中会明确警告, 不会静默降级。
// ============================================================================

include { FASTQC        } from '../../modules/local/qc/fastqc.nf'
include { FASTP         } from '../../modules/local/qc/fastp.nf'
include { FASTP_SUMMARY } from '../../modules/local/qc/fastp_summary.nf'
include { HOST_REMOVAL  } from '../../modules/local/qc/host_removal.nf'

workflow PREPROCESSING {

    take:
    ch_raw_reads    // channel: [ val(meta), path(reads) ]  reads=[R1,R2] or [R1]

    main:
    ch_versions     = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // ---------------------------------------------------------------------
    // 1. 原始 reads 质量评估
    //    与 fastp 并行 —— FastQC 不产出下游数据, 仅生成报告。
    // ---------------------------------------------------------------------
    if (!params.skip_fastqc) {
        FASTQC(ch_raw_reads)
        ch_versions      = ch_versions.mix(FASTQC.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { meta, zip -> zip })
    }

    // ---------------------------------------------------------------------
    // 2. 接头修剪与质量过滤
    // ---------------------------------------------------------------------
    FASTP(ch_raw_reads)
    ch_versions      = ch_versions.mix(FASTP.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map { meta, json -> json })

    // fastp 指标汇总表 (跨样本, 因此 collect)
    // 解析脚本以 file() 定位到本文件所在仓库, 与调用入口的 projectDir 无关
    ch_fastp_parser = file("${moduleDir}/../../bin/parse_fastp_json.py", checkIfExists: true)
    FASTP_SUMMARY(FASTP.out.json.map { meta, json -> json }.collect(), ch_fastp_parser)

    // ---------------------------------------------------------------------
    // 3. 宿主去除 (条件执行)
    // ---------------------------------------------------------------------
    def run_host_removal = params.host_index && !params.skip_host_removal

    if (run_host_removal) {
        // Bowtie2 索引由多个文件组成 (.1.bt2 ... .rev.2.bt2)。
        // 全部收集为一个 value channel, 供每个样本复用而不重复暂存。
        ch_host_index = Channel
            .fromPath("${params.host_index}*.{bt2,bt2l}", checkIfExists: true)
            .collect()

        HOST_REMOVAL(FASTP.out.reads, ch_host_index)

        ch_clean_reads   = HOST_REMOVAL.out.reads
        ch_versions      = ch_versions.mix(HOST_REMOVAL.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(HOST_REMOVAL.out.log.map { meta, log -> log })
    }
    else {
        // 明确告知跳过原因, 避免误以为已完成宿主去除
        if (params.skip_host_removal) {
            log.warn "宿主去除已跳过 (--skip_host_removal)。下游使用的是 fastp 清洁 reads, 未去除宿主序列。"
        }
        else {
            log.warn "未提供 --host_index, 宿主去除已跳过。下游使用的是 fastp 清洁 reads, 未去除宿主序列。"
        }

        ch_clean_reads = FASTP.out.reads
    }

    emit:
    // 下游唯一的 reads 来源 —— read-based 与 assembly 分支共用
    // [ val(meta), path(reads) ]  reads=[R1,R2] (PE) 或 [R1] (SE)
    reads         = ch_clean_reads

    // QC 产物 (供 MultiQC 与后续汇总)
    fastp_json    = FASTP.out.json            // [ val(meta), path(json) ]
    fastp_summary = FASTP_SUMMARY.out.tsv     // path(tsv)
    multiqc_files = ch_multiqc_files          // path(*)
    versions      = ch_versions               // path(versions.yml)
}
