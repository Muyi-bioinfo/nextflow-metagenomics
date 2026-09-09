// ============================================================================
// DIAMOND_BLASTP — 每个代表 MAG 的蛋白质组 blastp 比对 NR (Phase 12)
//
// 输入:  tuple(meta, mag_id, proteins.faa) —— GENE_PREDICTION.out.proteins
// 输出:  tuple(meta, mag_id, <mag_id>.diamond.tsv, proteins.faa) —— 原始
//        命中表 (faa 透传, 供 DIAMOND_SUMMARY 提取/校验 gene id)
//
// 数据库: params.diamond_db (NR 蛋白库 .dmnd), 缺失守卫在 ANNOTATION
//         子工作流。outfmt 6 固定 12 列, 解析契约见 bin/diamond_summary.py;
//         diamond_args 不要覆盖 --db / --query / --out / --outfmt。
// ============================================================================

process DIAMOND_BLASTP {
    tag "$mag_id"
    // Phase 17 资源审计: NR 库数十 GB, blastp 需 16GB+ 内存与多线程,
    // process_single (1 cpu/2GB/2h) 明显不足
    label 'process_high'

    conda "bioconda::diamond=2.1.11"
    container 'quay.io/biocontainers/diamond:2.1.11--h5ca1c30_2'

    publishDir { "${params.outdir}/${params.batch_id}/12_annotation/diamond" }, mode: 'copy', pattern: '*.diamond.tsv'

    input:
    tuple val(meta), val(mag_id), path(faa)

    output:
    tuple val(meta), val(mag_id), path("*.diamond.tsv"), path(faa), emit: hits
    path "versions.yml", emit: versions

    script:
    """
    diamond blastp \\
        --db ${params.diamond_db} \\
        --query ${faa} \\
        --out ${mag_id}.diamond.tsv \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \\
        --evalue 1e-5 \\
        --threads ${task.cpus} \\
        ${params.diamond_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | sed 's/^diamond version //')
    END_VERSIONS
    """

    stub:
    """
    # 占位输出仅验证通道形状与文件命名约定, 内容不代表真实比对结果。
    # 空命中表: PRODIGAL stub 的 faa 为空, DIAMOND_SUMMARY 走真实解析脚本,
    # raw 与 faa 的空/非空须一致 (否则未知 gene id 报错), 故不伪造命中行。
    touch ${mag_id}.diamond.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: 2.1.11
    END_VERSIONS
    """
}
