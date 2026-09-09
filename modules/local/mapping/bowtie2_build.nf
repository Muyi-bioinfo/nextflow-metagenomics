// ============================================================================
// BOWTIE2_BUILD — 为组装 contigs 构建 Bowtie2 索引
//
// 职责: 把一个组装单元的 contigs 变成可比对的参考索引。仅此一件事。
//
// 输入:  tuple(meta, contigs)   meta 来自 ASSEMBLY, 携带 assembler /
//                               assembly_mode / samples 字段
// 输出:  tuple(meta, bt2_index) 索引目录, 内含 <id>.<assembler>.*.bt2
//
// ─── 与 Phase 3 宿主索引的区别 ───────────────────────────────────────────
// Phase 3 的宿主索引是**外部输入** (--host_index 指向用户预先建好的索引);
// 这里的索引是**流程内产物** —— 参考序列就是本次组装出的 contigs, 每个组装
// 单元各建一份。因此没有任何路径需要硬编码: 参考序列从通道来。
//
// ─── 为什么输出成目录 ─────────────────────────────────────────────────────
// 一套 Bowtie2 索引是 6 个文件 (.1.bt2 ... .rev.2.bt2)。打包成单个目录后,
// 通道里只需传一个 path, 下游 process 暂存时也不会与 reads 等文件混在一起。
// 索引前缀由下游在运行时从 *.rev.1.bt2 反推 (见 bowtie2_map.nf), 不做约定式
// 拼接 —— 与 host_removal.nf 的做法保持一致。
//
// ─── 索引默认不发布 ───────────────────────────────────────────────────────
// 索引是可再生的中间产物, 体积与组装规模同量级。需要复用时用
// --save_bowtie2_index 显式打开。
// ============================================================================

process BOWTIE2_BUILD {
    tag "${meta.id}.${meta.assembler}"
    label 'process_medium'

    conda "bioconda::bowtie2=2.5.5"
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'

    publishDir { "${params.outdir}/${params.batch_id}/06_mapping/index/${meta.id}.${meta.assembler}" },
        mode: 'copy', pattern: "bt2_index", enabled: params.save_bowtie2_index

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("bt2_index"), emit: index
    path "versions.yml",                emit: versions

    script:
    def prefix = "${meta.id}.${meta.assembler}"
    """
    if [ ! -s ${contigs} ]; then
        echo "ERROR: contigs 文件为空 (${contigs}), 无法建立索引。" >&2
        echo "       上游组装应当在无 contig 时即报错, 请检查 05_assembly/ 的日志。" >&2
        exit 1
    fi

    mkdir -p bt2_index

    # contigs 为 gzip 压缩: bowtie2-build 原生支持 .gz 输入, 无需解压
    bowtie2-build \\
        --threads ${task.cpus} \\
        ${params.bowtie2_build_args} \\
        ${contigs} \\
        bt2_index/${prefix}

    # 以 .rev.1.bt2* 为锚点校验索引完整 —— 该后缀在一套索引中唯一
    # (小索引产出 .bt2, 大基因组产出 .bt2l, 二者存其一即可)
    if [ ! -e bt2_index/${prefix}.rev.1.bt2 ] && [ ! -e bt2_index/${prefix}.rev.1.bt2l ]; then
        echo "ERROR: bowtie2-build 未产出完整索引 (缺少 ${prefix}.rev.1.bt2/.bt2l)。" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$( bowtie2 --version 2>&1 | head -n 1 | sed 's/^.*bowtie2-align-s version //; s/ .*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}.${meta.assembler}"
    """
    mkdir -p bt2_index
    touch bt2_index/${prefix}.1.bt2 bt2_index/${prefix}.2.bt2 \\
          bt2_index/${prefix}.3.bt2 bt2_index/${prefix}.4.bt2 \\
          bt2_index/${prefix}.rev.1.bt2 bt2_index/${prefix}.rev.2.bt2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: 2.5.5
    END_VERSIONS
    """
}
