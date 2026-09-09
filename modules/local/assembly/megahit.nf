// ============================================================================
// MEGAHIT — 从头组装 (de novo assembly), 默认组装器
//
// 职责: 把清洁 reads 组装成 contigs。仅此一件事 —— 组装质量评估是 QUAST 的
//       职责 (见 quast.nf), 与本模块无关。
//
// 输入:  tuple(meta, reads)   reads 为列表: SE=[R1], PE=[R1,R2]
//                             meta 已由 subworkflows/local/assembly.nf 扩展,
//                             携带 assembler / assembly_mode / samples 字段
// 输出:  contigs              ${meta.id}.megahit.contigs.fa.gz
//        log                  MEGAHIT 运行日志 (k-mer 各轮次的 contig 统计)
//
// ─── 与 metaSPAdes 的关系 ─────────────────────────────────────────────────
// 二者是**可选的替代组装器**, 不是流水线的两个串联步骤。由 params.assembler
// 选择, "both" 时并行各跑一次 (各自独立的 assembly 单元)。见 assembly.nf。
//
// ─── 为什么默认选 MEGAHIT ─────────────────────────────────────────────────
// 内存占用远低于 metaSPAdes (succinct de Bruijn graph), 对复杂宏基因组样本
// 更容易跑完; 代价是 contig 连续性通常略逊于 metaSPAdes。
//
// ─── 输出压缩 ─────────────────────────────────────────────────────────────
// contigs 一律 gzip: 真实样本的 final.contigs.fa 可达数 GB。下游 QUAST /
// bowtie2-build 均原生支持 .gz 输入, 无需解压。
//
// ─── 空组装 ───────────────────────────────────────────────────────────────
// 组装不出任何 contig 时显式 exit 1 —— 不产出空文件冒充成功。MEGAHIT 自身在
// 极低深度下也可能以 "Too few vertices" 非零退出, 二者都会让任务可见地失败。
// ============================================================================

process MEGAHIT {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::megahit=1.2.9"
    container 'quay.io/biocontainers/megahit:1.2.9--haf24da9_8'

    publishDir { "${params.outdir}/${params.batch_id}/05_assembly/megahit/${meta.id}" },
        mode: 'copy', pattern: "*.megahit.{contigs.fa.gz,log}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.megahit.contigs.fa.gz"), emit: contigs
    tuple val(meta), path("${meta.id}.megahit.log"),           emit: log
    path "versions.yml",                                       emit: versions

    script:
    def prefix     = "${meta.id}.megahit"
    def input_arg  = meta.single_end ? "-r ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    // --k-list 与 --presets 互斥 (preset 本身就是一组 k); 互斥性在 assembly.nf 中前置校验
    def kmer_arg   = params.megahit_k_list  ? "--k-list ${params.megahit_k_list}"   : ''
    def preset_arg = params.megahit_preset  ? "--presets ${params.megahit_preset}"  : ''
    // 默认 2: 只出现一次的 (k+1)-mer 视为测序错误而滤除。极低深度数据需设为 1
    def count_arg  = params.megahit_min_count ? "--min-count ${params.megahit_min_count}" : ''
    // MEGAHIT 的 -m 接受 0-1 的比例或字节数整数。传字节数, 使其受 Nextflow 的
    // 资源分配约束, 而不是按整机内存的比例自行决定。
    def mem_bytes  = task.memory ? task.memory.toBytes() : 0
    def mem_arg    = mem_bytes > 0 ? "-m ${mem_bytes}" : ''
    """
    # MEGAHIT 拒绝写入已存在的输出目录, 因此固定用一个任务内的新目录
    megahit \\
        ${input_arg} \\
        -o megahit_out \\
        --out-prefix ${prefix} \\
        --num-cpu-threads ${task.cpus} \\
        --min-contig-len ${params.megahit_min_contig_len} \\
        ${mem_arg} \\
        ${kmer_arg} \\
        ${preset_arg} \\
        ${count_arg} \\
        ${params.megahit_args}

    if [ ! -s megahit_out/${prefix}.contigs.fa ]; then
        echo "ERROR: MEGAHIT 未产出任何 contig (megahit_out/${prefix}.contigs.fa 为空)。" >&2
        echo "       常见原因: 测序深度过低 / reads 过短 / 上游过滤过严。" >&2
        echo "       可尝试: --megahit_min_count 1, 用 --megahit_k_list 降低起始 k, 或增加测序量。" >&2
        exit 1
    fi

    gzip -c megahit_out/${prefix}.contigs.fa > ${prefix}.contigs.fa.gz

    # MEGAHIT 的日志文件名随版本而异: 1.2.9 写成 <out-prefix>.log, 更早的版本
    # 写成固定的 megahit_out/log。两者都接受, 都不存在才算异常。
    if [ -f megahit_out/${prefix}.log ]; then
        cp megahit_out/${prefix}.log ${prefix}.log
    elif [ -f megahit_out/log ]; then
        cp megahit_out/log ${prefix}.log
    else
        echo "ERROR: 未找到 MEGAHIT 运行日志 (megahit_out/${prefix}.log 或 megahit_out/log)。" >&2
        ls -1 megahit_out >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$( megahit --version 2>&1 | sed 's/^MEGAHIT v//' )
    END_VERSIONS
    """

    stub:
    """
    echo ">${meta.id}_contig_1" | gzip > ${meta.id}.megahit.contigs.fa.gz
    touch ${meta.id}.megahit.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: 1.2.9
    END_VERSIONS
    """
}
