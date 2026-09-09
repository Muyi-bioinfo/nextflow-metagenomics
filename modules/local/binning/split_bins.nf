// ============================================================================
// SPLIT_BINS — 将分箱器的 bin 目录拆分为带稳定 MAG ID 的独立 FASTA
//
// 职责: 把分箱器输出的 bin 目录 (内含 <prefix>.<N>.fa) 拆分成独立命名的
//       MAG FASTA, 并输出一张 per-MAG 统计表。
//
// 输入:  tuple(meta, bins_dir)
//          meta     组装单元的 meta (id / assembler / assembly_mode / samples)
//          bins_dir 分箱器输出目录 (METABAT2 / MaxBin2 / ...)
// 输出:  tuple(meta, mag_id, mag_fasta)   每个 MAG 一条 (flatMap 展开)
//        summary_tsv                       per-MAG 统计表 (n_contigs / total_bp / GC / ...)
//
// ─── MAG ID 的构成与稳定性 ────────────────────────────────────────────────
//     mag_id = <组装单元 id>.<assembler>.<binner>.<bin 序号, 三位补零>
//     例:      S01.megahit.metabat2.001
//
// 四个字段各自解决一类撞名:
//   组装单元 id   不同样本 (或 coassembly 组) 的 bin
//   assembler     --assembler both 时同一样本有两套独立 contigs
//   binner        V2 接入 MaxBin2 / CONCOCT / DAS Tool 后同一套 contigs 有多份分箱
//   bin 序号      同一次分箱内的不同 bin
//
// 序号**沿用分箱器自己给出的编号**, 不重新连续编号。详见 bin/split_bins.py。
//
// ─── 输出通道的形状 ───────────────────────────────────────────────────────
// 本 process 的 output 是 `tuple val(meta), path("mags/*.fa"), emit: mags`
// (bins 目录有 N 个 bin 时就是 N 个文件), 然后 emit 块用 flatMap 展开成
// `tuple(meta, mag_id, mag_fasta)` —— 每个 MAG 一条记录, 方便 Phase 8
// CheckM2 逐 MAG 处理。
//
// ─── 0 个 bin 时的处理 ────────────────────────────────────────────────────
// bin/split_bins.py 在 0 个 bin 时 exit 0, 写出仅含表头的 summary。本 process
// 产出一个空的 mags/ 目录 + 表头 summary —— 空通道自然不会触发下游任务,
// 不需要额外守卫。
// ============================================================================

process SPLIT_BINS {
    tag "${meta.id}.${meta.assembler}.${params.binner}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container 'quay.io/biocontainers/python:3.12'

    // pattern 保留目录结构会与路径末尾的 mags/ 叠加成 mags/mags/ 双层嵌套
    // (Phase 19 修复): saveAs 取 basename 平铺到 07_binning/mags/ 单层。
    // 仅改发布位置, 内部通道 emit (work/ 布局) 与下游 flatMap 不受影响。
    publishDir { "${params.outdir}/${params.batch_id}/07_binning/mags" },
        mode: 'copy', pattern: "mags/*.fa", saveAs: { it.split('/')[-1] }
    publishDir { "${params.outdir}/${params.batch_id}/07_binning" },
        mode: 'copy', pattern: "*_summary.tsv"

    input:
    tuple val(meta), path(bins_dir)
    path splitter   // bin/split_bins.py 作为显式输入, 使本 process 可从任意入口复用

    output:
    tuple val(meta), path("mags/*.fa"), optional: true, emit: mags   // 0 bin 时为空, flatMap 自然是空通道
    path "${meta.id}.${meta.assembler}.${params.binner}_summary.tsv", emit: summary
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}.${meta.assembler}.${params.binner}"
    """
    mkdir -p mags

    python3 ${splitter} \\
        --bins-dir ${bins_dir} \\
        --unit-id ${meta.id} \\
        --assembler ${meta.assembler} \\
        --assembly-mode ${meta.assembly_mode} \\
        --binner ${params.binner} \\
        --outdir mags \\
        --summary ${prefix}_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python3 --version | sed 's/Python //' )
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}.${meta.assembler}.${params.binner}"
    """
    mkdir -p mags
    echo ">${meta.id}.${meta.assembler}.${params.binner}.001_contig1" > mags/${meta.id}.${meta.assembler}.${params.binner}.001.fa
    echo "ACGTACGT" >> mags/${meta.id}.${meta.assembler}.${params.binner}.001.fa

    printf "mag_id\\tassembly_unit\\tassembler\\tassembly_mode\\tbinner\\tsource_bin\\tn_contigs\\ttotal_bp\\tlargest_contig_bp\\tmean_contig_bp\\tgc_percent\\tn_bases\\n" > ${prefix}_summary.tsv
    printf "${meta.id}.${meta.assembler}.${params.binner}.001\\t${meta.id}\\t${meta.assembler}\\t${meta.assembly_mode}\\t${params.binner}\\t${prefix}.1.fa\\t1\\t8\\t8\\t8.0\\t50.00\\t0\\n" >> ${prefix}_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
