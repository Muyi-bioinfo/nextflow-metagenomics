// ============================================================================
// CONTIG_DEPTH — contig 级覆盖度矩阵 (MetaBAT2 输入)
//
// 职责: 把一个组装单元的全部样本 BAM 汇总成 contig × 样本 的深度矩阵。
//
// 输入:  tuple(meta, bams, bais)
//          meta  组装单元的 meta (id / assembler / assembly_mode / samples)
//          bams  比对到**同一套 contigs** 的全部排序 BAM
// 输出:  tuple(meta, depth.txt)   MetaBAT2 --abdFile 直接可用
//
// ─── 为什么用 jgi_summarize_bam_contig_depths 而不是 samtools depth ───────
// MetaBAT2 需要的不是逐碱基深度, 而是每条 contig 的**平均深度与方差**, 且
// 会剔除 contig 两端 (默认各 75 bp) 的边缘效应、按比对质量加权。该工具随
// MetaBAT2 一同发布, 输出格式与 MetaBAT2 的解析器逐列对应 —— 自己用
// samtools depth 拼一份等价矩阵既易错也无必要。
//
// ─── 输出格式 (MetaBAT2 --abdFile 要求) ───────────────────────────────────
//   contigName  contigLen  totalAvgDepth  <bam1>  <bam1>-var  <bam2>  ...
// 即固定 3 列 + 每个 BAM 两列 (均值/方差)。script 末尾对此做实际校验:
// 列数与 BAM 数不匹配、或没有任何 contig 记录, 都直接失败, 不产出一个
// 「格式看起来对但内容是空的」文件流到 Phase 7。
//
// ─── 多样本列的意义 ───────────────────────────────────────────────────────
// 共变异 (co-abundance) 是分箱的核心信号之一: 同一个基因组的 contigs 在不同
// 样本间深度应当同步涨落。single 模式下每个组装单元只有本样本一个 BAM, 因此
// 只有单列深度; coassembly 接入后同一套 contigs 会有该组全部样本的 BAM,
// MetaBAT2 才能用上共变异信号 —— 这也是本模块按**组装单元**聚合 BAM
// (而非按样本逐个出文件) 的原因。
//
// ─── 默认不做长度/深度过滤 ────────────────────────────────────────────────
// --minContigLength / --minContigDepth 默认沿用工具本身的默认值, 保持深度
// 矩阵完整: contig 长度过滤是 MetaBAT2 自己的职责 (-m, 默认 2500), 在此处
// 提前砍掉会让同一份深度文件无法再被其他下游 (如 Phase 13 丰度计算) 复用。
// ============================================================================

process CONTIG_DEPTH {
    tag "${meta.id}.${meta.assembler}"
    label 'process_medium'

    conda "bioconda::metabat2=2.18"
    container 'quay.io/biocontainers/metabat2:2.18--h38e344b_2'

    publishDir { "${params.outdir}/${params.batch_id}/06_mapping/depth" },
        mode: 'copy', pattern: "*.depth.txt"

    input:
    tuple val(meta), path(bams), path(bais)

    output:
    tuple val(meta), path("${meta.id}.${meta.assembler}.depth.txt"), emit: depth
    path "versions.yml",                                             emit: versions

    script:
    def prefix        = "${meta.id}.${meta.assembler}"
    def bam_list      = bams instanceof List ? bams : [ bams ]
    def n_bam         = bam_list.size()
    def min_len_arg   = params.coverage_min_contig_len   ? "--minContigLength ${params.coverage_min_contig_len}" : ''
    def min_depth_arg = params.coverage_min_depth        ? "--minContigDepth ${params.coverage_min_depth}"       : ''
    """
    jgi_summarize_bam_contig_depths \\
        --outputDepth ${prefix}.depth.txt \\
        ${min_len_arg} \\
        ${min_depth_arg} \\
        ${bam_list.join(' ')}

    # ---- 实际校验输出确实是 MetaBAT2 可用的深度矩阵 ----
    if [ ! -s ${prefix}.depth.txt ]; then
        echo "ERROR: 未产出覆盖度文件 (${prefix}.depth.txt 为空)。" >&2
        exit 1
    fi

    N_COL=\$( head -n 1 ${prefix}.depth.txt | awk -F'\\t' '{print NF}' )
    N_ROW=\$( awk 'NR > 1' ${prefix}.depth.txt | wc -l )
    EXPECTED_COL=\$(( 3 + 2 * ${n_bam} ))

    if [ "\$N_COL" -ne "\$EXPECTED_COL" ]; then
        echo "ERROR: 覆盖度矩阵列数为 \$N_COL, 与 ${n_bam} 个 BAM 应有的 \$EXPECTED_COL 列不符。" >&2
        echo "       MetaBAT2 要求: contigName contigLen totalAvgDepth + 每个 BAM 两列 (均值/方差)。" >&2
        head -n 1 ${prefix}.depth.txt >&2
        exit 1
    fi

    if [ "\$N_ROW" -lt 1 ]; then
        echo "ERROR: 覆盖度矩阵没有任何 contig 记录 (只有表头)。" >&2
        echo "       常见原因: BAM 未坐标排序, 或没有任何 read 比对上 contigs。" >&2
        exit 1
    fi

    echo "覆盖度矩阵: \$N_ROW 条 contig × ${n_bam} 个样本 BAM (共 \$N_COL 列)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jgi_summarize_bam_contig_depths: \$( jgi_summarize_bam_contig_depths 2>&1 | head -n 1 | sed 's/^jgi_summarize_bam_contig_depths //; s/ .*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix   = "${meta.id}.${meta.assembler}"
    def bam_list = bams instanceof List ? bams : [ bams ]
    // 与真实输出同构: 固定 3 列 + 每个 BAM 两列 (均值/方差)
    def header   = ([ 'contigName', 'contigLen', 'totalAvgDepth' ] +
                    bam_list.collect { [ "${it}", "${it}-var" ] }.flatten()).join('\t')
    def row      = ([ "${meta.id}_contig_1", '1000', '10.00' ] +
                    bam_list.collect { [ '10.00', '1.00' ] }.flatten()).join('\t')
    """
    printf '%s\\n' '${header}' >  ${prefix}.depth.txt
    printf '%s\\n' '${row}'    >> ${prefix}.depth.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jgi_summarize_bam_contig_depths: 2.18
    END_VERSIONS
    """
}
