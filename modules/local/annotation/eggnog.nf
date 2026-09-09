// ============================================================================
// EGGNOG_MAPPER — 每个代表 MAG 的蛋白质组 eggNOG 功能注释 (Phase 12)
//
// 输入:  tuple(meta, mag_id, proteins.faa) —— GENE_PREDICTION.out.proteins
// 输出:  tuple(meta, mag_id, <mag_id>.emapper.annotations, proteins.faa)
//        (faa 透传, 供 EGGNOG_SUMMARY 提取/校验 gene id)
//
// 数据库: params.eggnog_db (eggNOG 数据目录, 含 eggnog.db 与
//         eggnog_proteins.dmnd), 缺失守卫在 ANNOTATION 子工作流。
//         GO 证据码经 --go_evidence 产出; 解析契约见 bin/eggnog_summary.py
//         (KO 多值逗号串/GO/COG 类别字母原样保留)。eggnog_args 不要覆盖
//         -i / -o / --output_dir / -m / --data_dir / --cpu。
// ============================================================================

process EGGNOG_MAPPER {
    tag "$mag_id"
    // Phase 17 资源审计: eggNOG 库 40+ GB, emapper 内部跑 diamond, 需
    // 16GB+ 内存与多线程 (--cpu), process_single 明显不足
    label 'process_high'

    conda "bioconda::eggnog-mapper=2.1.12"
    container 'quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0'

    publishDir { "${params.outdir}/${params.batch_id}/12_annotation/eggnog" }, mode: 'copy', pattern: '*.emapper.{annotations,seed_orthologs}'

    input:
    tuple val(meta), val(mag_id), path(faa)

    output:
    tuple val(meta), val(mag_id), path("*.emapper.annotations"), path(faa), emit: annotations
    path "*.emapper.seed_orthologs", emit: seed_orthologs, optional: true
    path "versions.yml", emit: versions

    script:
    """
    emapper.py \\
        -i ${faa} \\
        -o ${mag_id} \\
        --output_dir . \\
        -m diamond \\
        --data_dir ${params.eggnog_db} \\
        --go_evidence non-electronic \\
        --cpu ${task.cpus} \\
        ${params.eggnog_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$(emapper.py --version 2>&1 | tail -1)
    END_VERSIONS
    """

    stub:
    """
    # 占位输出仅验证通道形状与文件命名约定, 内容不代表真实注释结果。
    # 空 annotations: PRODIGAL stub 的 faa 为空, EGGNOG_SUMMARY 走真实解析
    # 脚本, raw 与 faa 的空/非空须一致, 故不伪造注释行。
    touch ${mag_id}.emapper.annotations ${mag_id}.emapper.seed_orthologs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: 2.1.12
    END_VERSIONS
    """
}
