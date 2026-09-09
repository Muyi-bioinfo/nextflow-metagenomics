// ============================================================================
// METASPADES — 从头组装 (de novo assembly), MEGAHIT 的替代组装器
//
// 职责: 把清洁 reads 组装成 contigs。与 MEGAHIT 完全平级 —— 由
//       params.assembler 二选一 (或 "both" 时并行各跑一次), 绝非串联的第二步。
//
// 输入:  tuple(meta, reads)   reads=[R1,R2]; metaSPAdes 不支持单端 (见下)
// 输出:  contigs              ${meta.id}.metaspades.contigs.fa.gz
//        scaffolds            ${meta.id}.metaspades.scaffolds.fa.gz
//        graph                assembly graph (GFA), 默认不发布
//        log                  spades.log
//
// ─── 单端限制 ─────────────────────────────────────────────────────────────
// `spades.py --meta` 要求恰好一个双端 library —— 这是 metaSPAdes 的硬性限制,
// 不是本模块的取舍。单端样本在 assembly.nf 中被分流并明确警告; 本模块保留一个
// 兜底分支, 使直接调用时也能给出可读的失败原因而非 SPAdes 的晦涩报错。
//
// ─── contigs 与 scaffolds ─────────────────────────────────────────────────
//   contigs    仅由 de Bruijn 图路径拼出的连续序列
//   scaffolds  在 contigs 之上, 用配对信息跨过 gap 连接, 中间以 N 填充
// 下游分箱 (MetaBAT2) 使用 contigs —— scaffold 中的 N gap 会干扰四核苷酸频率
// 与覆盖度计算。scaffolds 一并保留供人工审阅, 不进入主数据流。
//
// ─── 内存 ─────────────────────────────────────────────────────────────────
// metaSPAdes 内存占用显著高于 MEGAHIT。--memory 以 GB 为单位, 取自 Nextflow
// 分配给本 task 的内存, 而非让 SPAdes 按整机内存自行决定。
// ============================================================================

process METASPADES {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::spades=4.3.0"
    container 'quay.io/biocontainers/spades:4.3.0--hde4eca7_1'

    publishDir { "${params.outdir}/${params.batch_id}/05_assembly/metaspades/${meta.id}" },
        mode: 'copy', pattern: "*.metaspades.{contigs.fa.gz,scaffolds.fa.gz,log}"
    publishDir { "${params.outdir}/${params.batch_id}/05_assembly/metaspades/${meta.id}" },
        mode: 'copy', pattern: "*.metaspades.graph.gfa.gz", enabled: params.save_assembly_graph

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.metaspades.contigs.fa.gz"),   emit: contigs
    tuple val(meta), path("${meta.id}.metaspades.scaffolds.fa.gz"), emit: scaffolds, optional: true
    tuple val(meta), path("${meta.id}.metaspades.graph.gfa.gz"),    emit: graph,     optional: true
    tuple val(meta), path("${meta.id}.metaspades.log"),             emit: log
    path "versions.yml",                                            emit: versions

    script:
    def prefix   = "${meta.id}.metaspades"
    def kmer_arg = params.metaspades_k ? "-k ${params.metaspades_k}" : ''
    // --memory 单位为 GB 且必须为正整数; 分配不足 1 GB 时兜底为 1
    def mem_gb   = task.memory ? Math.max(1L, task.memory.toGiga()) : 1L

    if (meta.single_end) {
        """
        echo "ERROR: metaSPAdes (spades.py --meta) 仅支持双端 library, 而样本 ${meta.id} 为单端。" >&2
        echo "       请对单端样本使用 --assembler megahit。" >&2
        exit 1
        """
    } else {
        """
        spades.py \\
            --meta \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -o spades_out \\
            --threads ${task.cpus} \\
            --memory ${mem_gb} \\
            ${kmer_arg} \\
            ${params.metaspades_args}

        if [ ! -s spades_out/contigs.fasta ]; then
            echo "ERROR: metaSPAdes 未产出任何 contig (spades_out/contigs.fasta 为空)。" >&2
            echo "       常见原因: 测序深度过低 / reads 过短 / 上游过滤过严。" >&2
            echo "       可尝试: 用 --metaspades_k 指定更小的 k 值列表, 或增加测序量。" >&2
            exit 1
        fi

        gzip -c spades_out/contigs.fasta > ${prefix}.contigs.fa.gz
        cp spades_out/spades.log ${prefix}.log

        # scaffolds 与 assembly graph 并非所有运行都产出, 存在才转存
        if [ -s spades_out/scaffolds.fasta ]; then
            gzip -c spades_out/scaffolds.fasta > ${prefix}.scaffolds.fa.gz
        fi
        if [ -s spades_out/assembly_graph_with_scaffolds.gfa ]; then
            gzip -c spades_out/assembly_graph_with_scaffolds.gfa > ${prefix}.graph.gfa.gz
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            spades: \$( spades.py --version 2>&1 | sed 's/^.*v//' )
        END_VERSIONS
        """
    }

    stub:
    """
    echo ">${meta.id}_contig_1" | gzip > ${meta.id}.metaspades.contigs.fa.gz
    echo ">${meta.id}_scaffold_1" | gzip > ${meta.id}.metaspades.scaffolds.fa.gz
    touch ${meta.id}.metaspades.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: 4.3.0
    END_VERSIONS
    """
}
