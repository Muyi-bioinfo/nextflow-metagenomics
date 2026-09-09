// ============================================================================
// PRODIGAL — 基因预测 (Phase 11)
//
// 输入:  tuple(meta, mag_id, mag_fasta) —— 每个代表 MAG 各调用一次 (来自
//        DREP.out.representatives; --skip_dereplication 时代表集 = 全部
//        qualified MAG, 该通道恒等传递, 永不空)
// 输出:  tuple(meta, mag_id, genes.fna)   / proteins.faa / gff
//
// -p meta 宏基因组模式; 无外部数据库依赖 (独立二进制)。
// 文件命名以 mag_id 为准 —— mag_id 本身 = <组装单元 id>.<assembler>.<binner>.
// <bin 序号> (见 bin/split_bins.py), 组装单元 id 即 meta.id, 故文件名天然与
// meta.id + mag_id 绑定 (CheckM2 同款约定)。GFF / FASTA 内部的基因 ID 由
// Prodigal 从原 MAG contig 头派生 (${seq_header}_N), 可回溯到 MAG 内的 contig。
// ============================================================================

process PRODIGAL {
    tag "${meta.id}.${mag_id}"
    label 'process_single'

    conda "bioconda::prodigal=2.6.3"
    container 'quay.io/biocontainers/prodigal:2.6.3--h577a1d6_11'

    publishDir { "${params.outdir}/${params.batch_id}/11_gene_prediction" },
        mode: 'copy', pattern: '*.{fna,faa,gff}'

    input:
    tuple val(meta), val(mag_id), path(mag_fasta)

    output:
    tuple val(meta), val(mag_id), path("${mag_id}.genes.fna"),    emit: genes
    tuple val(meta), val(mag_id), path("${mag_id}.proteins.faa"), emit: proteins
    tuple val(meta), val(mag_id), path("${mag_id}.gff"),          emit: gff
    path "versions.yml", emit: versions

    script:
    """
    prodigal \\
        -i ${mag_fasta} \\
        -a ${mag_id}.proteins.faa \\
        -d ${mag_id}.genes.fna \\
        -f gff \\
        -o ${mag_id}.gff \\
        -p meta \\
        ${params.prodigal_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prodigal: 2.6.3
    END_VERSIONS
    """

    stub:
    """
    # 占位输出仅验证通道形状与文件命名约定, 内容不代表真实预测结果
    touch ${mag_id}.genes.fna ${mag_id}.proteins.faa ${mag_id}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prodigal: 2.6.3
    END_VERSIONS
    """
}
