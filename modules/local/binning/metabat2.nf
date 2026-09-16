// ============================================================================
// METABAT2 — MAG 分箱 (Phase 7)
//
// 职责: 把一个组装单元的 contigs 分箱 (binning), 产出 bin 目录。
//
// 输入:  tuple(meta, contigs, depth)
//          meta    组装单元的 meta (id / assembler / assembly_mode / samples)
//          contigs 组装器产出的 FASTA (gzip 或非 gzip 均可)
//          depth   Phase 6 的覆盖度矩阵 (jgi_summarize_bam_contig_depths 输出)
// 输出:  tuple(meta, bins_dir)   bins 目录 (内含 <prefix>.<N>.fa)
//        log                      MetaBAT2 运行日志
//
// ─── MetaBAT2 的运行模式 ──────────────────────────────────────────────────
// MetaBAT2 以**组装单元级**运行 —— 一套 contigs + 该单元的全样本深度矩阵 →
// 一次分箱。不是逐样本分箱: 共变异信号 (同一基因组的 contigs 在不同样本间
// 深度同步涨落) 是 MetaBAT2 的核心依据之一, single 模式下只有单列深度因此
// 只能靠组成 (tetranucleotide frequency); coassembly 下同一套 contigs 有
// 该组全部样本的深度, MetaBAT2 才能用上共变异。
//
// ─── 输出的是「bin 目录」而非「MAG 文件列表」 ─────────────────────────
// MetaBAT2 产出一个目录, 内含 <prefix>.1.fa / <prefix>.2.fa / ..., 可能还有
// <prefix>.unbinned.fa / <prefix>.tooShort.fa / <prefix>.lowDepth.fa 等
// 非 bin 产物。这些文件名只在单个组装单元内唯一 —— 一旦多个样本 / 多个
// 组装器 / 多个分箱器的结果汇到一起 (Phase 8 CheckM2 / Phase 9 dRep),
// 就会互相撞名。
//
// 本 process 只负责「运行 MetaBAT2, 把 bin 目录作为整体输出」; 拆分成带
// 全局唯一 MAG ID 的独立 FASTA 是下游 SPLIT_BINS 的职责 (见 split_bins.nf)。
//
// ─── 0 个 bin 是合法结果 ──────────────────────────────────────────────────
// 低复杂度或低深度数据下 MetaBAT2 形成不了满足最小 bin 尺寸 (默认 200 kb)
// 的簇是真实结果, 不是失败。本 process 正常退出 (exit 0), 只在 stderr 说明。
// 虚构一个空 bin 才是错的。
//
// ─── --seed 与 MAG ID 稳定性 ──────────────────────────────────────────────
// MetaBAT2 内部使用随机采样 (bootstrap), 不固定种子则同样的输入可能产出
// 不同的 bin 划分 —— MAG ID 也就跟着变 (见 bin/split_bins.py 的文件头)。
// params.metabat2_seed 默认值应当是一个固定整数 (如 42), 而非 null / -1。
// ============================================================================

process METABAT2 {
    tag "${meta.id}.${meta.assembler}"
    label 'process_high'

    conda "bioconda::metabat2=2.18"
    container 'quay.io/biocontainers/metabat2:2.18--h38e344b_2'

    publishDir { "${params.outdir}/${params.batch_id}/07_binning/metabat2/${meta.id}.${meta.assembler}" },
        mode: 'copy', pattern: "bins"
    publishDir { "${params.outdir}/${params.batch_id}/07_binning/logs" },
        mode: 'copy', pattern: "*.log"

    input:
    tuple val(meta), path(contigs), path(depth)

    output:
    tuple val(meta), path("bins"), emit: bins
    path "${meta.id}.${meta.assembler}.metabat2.log", emit: log
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}.${meta.assembler}"
    // MetaBAT2 -m: 进入分箱的最小 contig 长度 (默认 2500); -s: 最小 bin 尺寸 (默认 200000)
    def min_contig_arg = params.metabat2_min_contig_len ? "-m ${params.metabat2_min_contig_len}" : ''
    def min_bin_arg    = params.metabat2_min_bin_size   ? "-s ${params.metabat2_min_bin_size}"   : ''
    def seed_arg       = params.metabat2_seed != null   ? "--seed ${params.metabat2_seed}"       : ''
    // --saveCls 输出每条 contig 的 bin 分配 (bin_id 或 0=unbinned), 可选; --saveProb
    // 输出每条 contig 的分配置信度; --noBinOut 只输出分配表不写 bin FASTA (调试用)。
    // 默认只要 bin FASTA, 不加这些开关。
    """
    # MetaBAT2 拒绝写入已存在的输出目录
    mkdir -p bins

    metabat2 \\
        -i ${contigs} \\
        -a ${depth} \\
        -o bins/${prefix} \\
        -t ${task.cpus} \\
        ${min_contig_arg} \\
        ${min_bin_arg} \\
        ${seed_arg} \\
        ${params.metabat2_args} \\
        > ${prefix}.metabat2.log 2>&1

    # ---- 0 个 bin 是合法结果: 不报错, 只说明 ----
    N_BINS=\$( ls -1 bins/${prefix}.[0-9]*.fa 2>/dev/null | wc -l )
    if [ "\$N_BINS" -eq 0 ]; then
        echo "未形成任何 bin (bins/ 目录为空, 或只有非 bin 产物如 .unbinned.fa)。" >&2
        echo "常见原因: 组装总长 < MetaBAT2 的最小 bin 尺寸 (默认 200 kb)," >&2
        echo "或最小 contig 长度阈值 (-m, 默认 2500) 过滤掉了全部 contigs。" >&2
        echo "低复杂度/低深度样本下这是真实结果。" >&2
    else
        echo "MetaBAT2 形成 \$N_BINS 个 bin -> bins/${prefix}.*.fa" >&2
        ls -1h bins/${prefix}.[0-9]*.fa >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabat2: \$( metabat2 --help 2>&1 | grep -i 'version' | sed 's/.*version //; s/ .*//' || echo "2.18" )
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}.${meta.assembler}"
    """
    mkdir -p bins
    echo ">${meta.id}_bin1_contig1" > bins/${prefix}.1.fa
    echo "ACGTACGT" >> bins/${prefix}.1.fa
    touch ${prefix}.metabat2.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabat2: 2.18
    END_VERSIONS
    """
}
