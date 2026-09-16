// ============================================================================
// RGI 分支 (Phase 12) — CARD 抗性基因注释, 两个 process:
//
// RGI_LOAD  一次性把 CARD 数据库 (params.card_db, card.json) 经
//           `rgi load --card_json ... --local` 载入本地库。rgi main 本身
//           不接收数据库路径参数, 只读本地库 —— 若每个 MAG 各载一次,
//           并行任务会并发写同一索引 (损坏风险), 故载入单独成 process,
//           以 .done 标记建立 RGI_MAIN 的依赖。conda 环境下本地库目录
//           可写; 容器部署时可在镜像内预载入, card_db 仍由守卫强制提供。
// RGI_MAIN  每个代表 MAG 的蛋白质组 RGI 分析 (蛋白模式)。蛋白模式下
//           ORF_ID 即输入 proteins.faa 的序列头 = Prodigal gene id,
//           是 Phase 14 join 键一致的前提 (rgi main 的 contig 模式会
//           自预测 ORF, id 与 Prodigal 不一致, 故不用)。
//
// 输入:  tuple(meta, mag_id, proteins.faa) —— GENE_PREDICTION.out.proteins
// 输出:  tuple(meta, mag_id, <mag_id>.rgi.json, proteins.faa)
//
// rgi_args 不要覆盖 -i / -o / -t / -a / --local / --clean / --threads。
// ============================================================================

process RGI_LOAD {
    label 'process_single'

    conda "bioconda::rgi=6.0.8"
    container 'quay.io/biocontainers/rgi:6.0.8--pyh05cac1d_0'

    input:
    // 上游 proteins 通道的首个元素仅作调度门: 通道为空时 (如
    // --skip_gene_prediction) 不触发无谓的数据库载入。
    tuple val(meta), val(mag_id), path(faa)

    output:
    path "rgi_load.done", emit: db
    path "versions.yml", emit: versions

    script:
    """
    rgi load \\
        --card_json ${params.card_db} \\
        --local

    touch rgi_load.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(rgi main --version 2>&1 | tail -1)
    END_VERSIONS
    """

    stub:
    """
    # 占位: 不做真实载入 (stub-run 中 card_db 通常是占位路径), 仅产出
    # 依赖标记。
    touch rgi_load.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: 6.0.8
    END_VERSIONS
    """
}

process RGI_MAIN {
    tag "$mag_id"
    // Phase 17 资源审计: rgi main 蛋白模式跑 diamond 比对, 2GB/1cpu/2h 过紧,
    // 提至 process_medium (RGI_LOAD 仅解析 card.json 一次, 保持 process_single)
    label 'process_medium'

    conda "bioconda::rgi=6.0.8"
    container 'quay.io/biocontainers/rgi:6.0.8--pyh05cac1d_0'

    publishDir { "${params.outdir}/${params.batch_id}/12_annotation/rgi" }, mode: 'copy', pattern: '*.rgi.json'

    input:
    tuple val(meta), val(mag_id), path(faa), path(rgi_db_ready)

    output:
    tuple val(meta), val(mag_id), path("*.rgi.json"), path(faa), emit: annotations
    path "versions.yml", emit: versions

    script:
    """
    rgi main \\
        --input_sequence ${faa} \\
        --input_type protein \\
        --alignment_tool DIAMOND \\
        --output_file ${mag_id}.rgi \\
        --local \\
        --clean \\
        --threads ${task.cpus} \\
        ${params.rgi_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(rgi main --version 2>&1 | tail -1)
    END_VERSIONS
    """

    stub:
    """
    # 占位输出仅验证通道形状与文件命名约定, 内容不代表真实注释结果。
    # 空 JSON (0 字节): PRODIGAL stub 的 faa 为空, RGI_SUMMARY 走真实解析
    # 脚本, raw 与 faa 的空/非空须一致, 故不虚构命中条目。
    touch ${mag_id}.rgi.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: 6.0.8
    END_VERSIONS
    """
}
