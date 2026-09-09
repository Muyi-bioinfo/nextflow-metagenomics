// ============================================================================
// nextflow-metagenomics — Nextflow 入口 (main.nf)
//
// Phase 1: Repository Skeleton (薄入口)
//
// main.nf 仅负责:
//   1. 启用 DSL2
//   2. 调用 workflows/mag.nf 中的主 workflow
//
// 所有流程逻辑位于 workflows/mag.nf，模块位于 modules/{nf-core,local}/，
// 子工作流位于 subworkflows/{nf-core,local}/。
// ============================================================================

nextflow.enable.dsl = 2

include { mag } from './workflows/mag.nf'

workflow {

    // 调用主 workflow (后续分析流程在 mag 中按阶段接入)
    mag()
}
