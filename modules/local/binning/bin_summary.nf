// ============================================================================
// BIN_SUMMARY — 汇总各组装单元的 per-MAG 统计表
//
// 职责: 把各个组装单元的 <unit>.<assembler>.<binner>_summary.tsv 纵向拼接成
//       一张 bin_summary.tsv, 类似 assembly.nf 的 ASSEMBLY_SUMMARY。
//
// 输入:  path(summaries)   各组装单元的 summary TSV 列表 (collect 聚合)
// 输出:  bin_summary.tsv   全局汇总表 (去重表头, 保留全部数据行)
//
// ─── 为什么需要这个汇总 ───────────────────────────────────────────────────
// Phase 14 (整合结果) 需要一张全局的 MAG 基础统计表 (contig 数 / 总长 / GC)
// 作为后续拼接的基底 —— CheckM2 评分 / GTDB-Tk 分类 / CoverM 丰度都要按
// mag_id 与这张表 join。不汇总的话 Phase 14 要自己遍历 07_binning/ 下
// 散落的各单元 summary, 不如在 Phase 7 本地就做了。
//
// ─── 0 个 MAG 时的处理 ───────────────────────────────────────────────────
// 所有单元都没形成 bin 时, collect() 拿到的是 N 个仅含表头的 TSV。拼接后
// 输出也是仅含表头的 bin_summary.tsv, 不是错误 (与 bin/split_bins.py 的
// "0 个 bin 是合法结果" 语义一致)。
// ============================================================================

process BIN_SUMMARY {
    label 'process_single'

    conda "conda-forge::python=3.12"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/07_binning" },
        mode: 'copy'

    input:
    path summaries   // 各单元 summary TSV 的列表 (可能是空列表, 或只有表头)

    output:
    path "bin_summary.tsv", emit: tsv
    path "versions.yml", emit: versions

    script:
    def summary_list = summaries instanceof List ? summaries.collect { s -> "'${s}'" }.join(', ') : "'${summaries}'"
    """
    #!/usr/bin/env python3
    import sys
    from pathlib import Path

    summaries = sorted([Path(p) for p in [${summary_list}]])

    if not summaries:
        # --skip_binning 或整条通道为空 (极端情况): 写一个占位表
        with open("bin_summary.tsv", "w") as out:
            out.write("mag_id\\tassembly_unit\\tassembler\\tassembly_mode\\tbinner\\t"
                      "source_bin\\tn_contigs\\ttotal_bp\\tlargest_contig_bp\\t"
                      "mean_contig_bp\\tgc_percent\\tn_bases\\n")
        print("未找到任何 bin summary (可能 --skip_binning, 或所有单元都未形成 bin)。", file=sys.stderr)
    else:
        header_written = False
        n_mags = 0
        with open("bin_summary.tsv", "w") as out:
            for tsv in summaries:
                with tsv.open("r") as fh:
                    lines = fh.readlines()
                    if not lines:
                        continue
                    # 第一个文件: 写表头 + 数据行; 后续文件: 只写数据行
                    if not header_written:
                        out.write(lines[0])
                        header_written = True
                    for line in lines[1:]:
                        out.write(line)
                        n_mags += 1

        if n_mags == 0:
            print("所有单元的 summary 都只有表头 (0 个 MAG)。", file=sys.stderr)
        else:
            print(f"汇总 {len(summaries)} 个单元, 共 {n_mags} 个 MAG -> bin_summary.tsv", file=sys.stderr)

    with open("versions.yml", "w") as v:
        v.write('"${task.process}":\\n')
        v.write('    python: "3.12"\\n')
    """

    stub:
    """
    printf "mag_id\\tassembly_unit\\tassembler\\tassembly_mode\\tbinner\\tsource_bin\\tn_contigs\\ttotal_bp\\tlargest_contig_bp\\tmean_contig_bp\\tgc_percent\\tn_bases\\n" > bin_summary.tsv
    printf "S01.megahit.metabat2.001\\tS01\\tmegahit\\tsingle\\tmetabat2\\tS01.megahit.1.fa\\t10\\t50000\\t8000\\t5000.0\\t45.2\\t0\\n" >> bin_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
