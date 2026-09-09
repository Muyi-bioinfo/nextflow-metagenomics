// ============================================================================
// KRAKEN2 — 基于 reads 的物种分类 (classification)
//
// 职责: 仅做 classification —— 为每条 read 指派分类单元, 并汇总为 report。
//       不做丰度估计, 那是 BRACKEN 的职责 (见 bracken.nf)。
//
// 输入:  tuple(meta, reads)   reads 为列表: SE=[R1], PE=[R1,R2]
//        path(db)             Kraken2 数据库目录 (含 hash.k2d/opts.k2d/taxo.k2d)
// 输出:  report               标准 Kraken2 报告 —— Bracken 的输入, MultiQC 可解析
//        classified_output    逐 read 的分类结果 (每行一条 read)
//
// ─── report 与 classification output 的区别 ────────────────────────────────
//   report            按分类树聚合的计数表 (百分比 / clade reads / taxon reads /
//                     rank / taxid / name), 每个分类单元一行 —— Bracken 只读这份。
//   classification    逐 read 一行 (C/U, read id, taxid, 长度, k-mer 命中映射),
//                     用于追溯单条 read 的判定依据。
//
// ─── 体积说明 ─────────────────────────────────────────────────────────────
// classification output 每条 read 一行, 真实样本可达数 GB。因此一律压缩, 且默认
// 不发布 (params.save_kraken2_output = false) —— 与 save_host_removed 同理:
// 通道中始终存在, 需要时开启发布即可, 不影响下游连接。
//
// ─── 数据库不复制 ─────────────────────────────────────────────────────────
// path(db) 由 Nextflow 以符号链接暂存而非复制 —— Kraken2 数据库动辄数十至数百 GB。
// ============================================================================

process KRAKEN2 {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::kraken2=2.17.1"
    container 'quay.io/biocontainers/kraken2:2.17.1--pl5321h077b44d_0'

    publishDir { "${params.outdir}/${params.batch_id}/03_taxonomy/kraken2" },
        mode: 'copy', pattern: "*.kraken2.{report.txt,log}"
    publishDir { "${params.outdir}/${params.batch_id}/03_taxonomy/kraken2/classifications" },
        mode: 'copy', pattern: "*.kraken2.output.txt.gz", enabled: params.save_kraken2_output

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("${meta.id}.kraken2.report.txt"),    emit: report
    tuple val(meta), path("${meta.id}.kraken2.output.txt.gz"), emit: classified_output
    tuple val(meta), path("${meta.id}.kraken2.log"),           emit: log
    path "versions.yml",                                       emit: versions

    script:
    def paired_arg  = meta.single_end ? '' : '--paired'
    def input_files = meta.single_end ? "${reads[0]}" : "${reads[0]} ${reads[1]}"
    def mmap_arg    = params.kraken2_memory_mapping ? '--memory-mapping' : ''
    """
    # 数据库完整性预检: 缺文件时 Kraken2 自身的报错不够直白, 这里先给出明确指引
    for f in hash.k2d opts.k2d taxo.k2d; do
        if [ ! -e "${db}/\$f" ]; then
            echo "ERROR: Kraken2 数据库不完整 —— 未找到 ${db}/\$f" >&2
            echo "       --kraken2_db 应指向包含 hash.k2d / opts.k2d / taxo.k2d 的目录。" >&2
            exit 1
        fi
    done

    kraken2 \\
        --db ${db} \\
        --threads ${task.cpus} \\
        --report ${meta.id}.kraken2.report.txt \\
        --output ${meta.id}.kraken2.output.txt \\
        --confidence ${params.kraken2_confidence} \\
        --minimum-base-quality ${params.kraken2_min_base_quality} \\
        --minimum-hit-groups ${params.kraken2_min_hit_groups} \\
        --gzip-compressed \\
        ${paired_arg} \\
        ${mmap_arg} \\
        ${input_files} \\
        2> >( tee ${meta.id}.kraken2.log >&2 )

    # 逐 read 结果体积大, 一律压缩后再交给下游/发布
    gzip -n ${meta.id}.kraken2.output.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$( kraken2 --version | head -n 1 | sed 's/^Kraken version //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.kraken2.report.txt ${meta.id}.kraken2.log
    echo | gzip > ${meta.id}.kraken2.output.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: 2.17.1
    END_VERSIONS
    """
}
