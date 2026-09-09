// ============================================================================
// SAMTOOLS_SORT — BAM 坐标排序 + 索引 + 比对统计
//
// 职责: 把 BOWTIE2_MAP 的未排序 BAM 变成「可被随机访问」的排序 BAM, 并产出
//       标准比对统计。这是覆盖度计算的**硬前提** ——
//       jgi_summarize_bam_contig_depths 按 contig 顺序流式扫描, 只接受坐标
//       排序的 BAM; 未排序输入会直接报错或给出错误深度。
//
// 输入:  tuple(meta, bam)      未排序 BAM
// 输出:  tuple(meta, bam, bai) 坐标排序 BAM + 索引
//        flagstat / idxstats / stats  (MultiQC 可解析)
//
// ─── 三个统计各自回答什么 ─────────────────────────────────────────────────
//   flagstat  总体比对率 —— reads 中有多大比例落回了自己的组装
//             (比例过低意味着组装未能代表样本主体, 分箱结果会失真)
//   idxstats  每条 contig 的比对 read 数 —— 快速定位异常高/低覆盖 contig
//   stats     插入片段分布、错配率等细节, 供 MultiQC 出图
//
// 这些是**比对质控**, 与 contig 覆盖度矩阵 (coverage.nf) 不同: 前者判断比对
// 本身是否可信, 后者是分箱算法的输入数据。
//
// ─── BAM 默认不发布 ───────────────────────────────────────────────────────
// 真实样本单个 BAM 可达数十 GB, 且下游 (覆盖度 / Phase 13 CoverM) 都从通道
// 直接取用, 无需落盘到结果目录。需要时用 --save_bam 打开。
// ============================================================================

process SAMTOOLS_SORT {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::samtools=1.24"
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'

    publishDir { "${params.outdir}/${params.batch_id}/06_mapping/bam" },
        mode: 'copy', pattern: "*.sorted.bam*", enabled: params.save_bam
    publishDir { "${params.outdir}/${params.batch_id}/06_mapping/logs" },
        mode: 'copy', pattern: "*.{flagstat,idxstats,stats}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.flagstat"),                                      emit: flagstat
    tuple val(meta), path("${meta.id}.idxstats"),                                      emit: idxstats
    tuple val(meta), path("${meta.id}.stats"),                                         emit: stats
    path "versions.yml",                                                               emit: versions

    script:
    // samtools sort 的 -m 是**每线程**上限, 且 -@ N 表示额外 N 个线程,
    // 峰值内存约 (N+1) × m。留 30% 余量, 避免刚好压在 task.memory 上被 OOM。
    def avail_mb   = task.memory ? task.memory.toMega() : 0
    def per_thread = avail_mb > 0 ? Math.max(768L, (avail_mb * 0.7 / (task.cpus + 1)) as long) : 768L
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -m ${per_thread}M \\
        -o ${meta.id}.sorted.bam \\
        ${bam}

    samtools index -@ ${task.cpus} ${meta.id}.sorted.bam

    samtools flagstat -@ ${task.cpus} ${meta.id}.sorted.bam > ${meta.id}.flagstat
    samtools idxstats ${meta.id}.sorted.bam                 > ${meta.id}.idxstats
    samtools stats -@ ${task.cpus} ${meta.id}.sorted.bam    > ${meta.id}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version | head -n 1 | sed 's/samtools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.bam ${meta.id}.sorted.bam.bai
    touch ${meta.id}.flagstat ${meta.id}.idxstats ${meta.id}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.24
    END_VERSIONS
    """
}
