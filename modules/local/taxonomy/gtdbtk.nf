// ============================================================================
// GTDBTK — MAG 物种分类 (Phase 10)
//
// GTDB-Tk classify_wf 是**批处理**工具: 本 process 把本阶段全部输入 MAG 聚合
// 为一个基因组目录, 一次调用完成分类 (逐 MAG 调用会在建树/比对步骤上
// 白白重复开销)。
//
// 输入:  genome_manifest.tsv (mag_id \t sample, 由 collectFile 物化)
//        + 输入 MAG FASTA 路径列表 (toSortedList 聚合) —— 执行顺序按最佳
//        实践为先去冗余后分类 (Phase 9 dRep), 故正常路径下这里是**代表 MAG
//        集合**; dRep 未接入时 (当前) 则为全部 qualified MAG。
// 输出:  gtdbtk.bac120.summary.tsv / gtdbtk.ar53.summary.tsv
//        (两者**恒存在**, 缺其一只表示该域无基因组, 为下游提供稳定接口)
//
// 映射回 (meta.id, mag_id) 的方式: 每个 MAG 以其 FASTA 文件名 (即 <mag_id>.fa,
// 由 Phase 7 的稳定 MAG ID 命名保证) 软链入 genome_dir, summary 的 user_genome
// 列即此文件名, taxonomy_summary 剥掉 .fa 扩展名后经 manifest 找回样本 ——
// MAG 身份在聚合之后不丢失。
//
// 注意: 不用 val 通道传 tuple 列表 —— Nextflow 的 collect 聚合会把 tuple
// 扁平化为 ArrayBag, process 内无法按三元组还原 (实测踩坑, 见 STATUS.md)。
// ============================================================================

process GTDBTK {
    tag "batch: ${params.batch_id}"
    label 'process_high'

    conda "bioconda::gtdbtk=2.7.2"
    container 'quay.io/biocontainers/gtdbtk:2.7.2--pyhdfd78af_1'

    // 汇总表平铺发布; 注意 Phase 7 的教训 —— publishDir 路径末尾与 pattern
    // 都带目录名会产生 mags/mags/ 式双层嵌套, 这里 pattern 不含目录前缀。
    publishDir { "${params.outdir}/${params.batch_id}/10_mag_taxonomy/gtdbtk" },
        mode: 'copy', pattern: "gtdbtk.*.summary.tsv"
    publishDir { "${params.outdir}/${params.batch_id}/10_mag_taxonomy" },
        mode: 'copy', pattern: "gtdbtk_out"

    input:
    path(manifest)
    path(fastas)

    output:
    path "gtdbtk.bac120.summary.tsv", emit: bac_summary
    path "gtdbtk.ar53.summary.tsv",  emit: ar_summary
    path "gtdbtk_out", type: 'dir', emit: raw
    path "versions.yml",             emit: versions

    script:
    def fasta_list = fastas instanceof List ? fastas : [fastas]
    def staging = fasta_list.collect { f -> "ln -s ${f} genome_dir/${f.getFileName()}" }.join('\n')
    """
    mkdir -p genome_dir

    # MAG ID 必须唯一 —— GTDB-Tk 以 genome 名为标识, 重名会让映射表失效
    if [ "\$(awk -F'\\t' 'NF>=2 {print \$1}' genome_manifest.tsv | sort | uniq -d | wc -l)" -ne 0 ]; then
        echo "ERROR: duplicate MAG ID in genome_manifest.tsv — GTDB-Tk genome names must be unique" >&2
        exit 1
    fi

    # 每个 MAG 以其 FASTA 文件名 (即 <mag_id>.fa) 软链入同一基因组目录
    ${staging}

    if [ -s genome_manifest.tsv ]; then
        export GTDBTK_DATA_PATH="${params.gtdbtk_db}"

        # 重试安全: classify_wf 不支持覆盖已有 out_dir
        rm -rf gtdbtk_out
        gtdbtk classify_wf \\
            --genome_dir genome_dir \\
            --out_dir gtdbtk_out \\
            -x fa \\
            --cpus ${task.cpus} \\
            ${params.gtdbtk_args}
    else
        echo "WARN: genome_manifest.tsv is empty — GTDB-Tk classification skipped" >&2
        mkdir -p gtdbtk_out
    fi

    # 细菌 (bac120) 与古菌 (ar53) summary 可能只产其一 —— 两者都保证在 task
    # 根目录存在 (缺其一即为空文件), 下游 TAXONOMY_SUMMARY 无需 optional 输入。
    cp gtdbtk_out/gtdbtk.bac120.summary.tsv gtdbtk.bac120.summary.tsv 2>/dev/null || : > gtdbtk.bac120.summary.tsv
    cp gtdbtk_out/gtdbtk.ar53.summary.tsv  gtdbtk.ar53.summary.tsv  2>/dev/null || : > gtdbtk.ar53.summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: 2.7.2
    END_VERSIONS
    """

    stub:
    """
    if [ -s genome_manifest.tsv ]; then
        # stub: 不真实调用 gtdbtk —— 由 manifest 生成假的 bac120 summary
        # (每 MAG 一行, 分类串含缺级空 token, 覆盖 -1 split 场景);
        # ar53 留空, 覆盖"两个 summary 缺其一"场景。
        mkdir -p gtdbtk_out
        printf 'user_genome\\taccession\\tclassification\\n' > gtdbtk_out/gtdbtk.bac120.summary.tsv
        while read -r mag_id sample; do
            [ -n "\$mag_id" ] || continue
            printf '%s.fa\\tGCA_stub_%s\\td__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__;g__;s__\\n' "\$mag_id" "\$mag_id" >> gtdbtk_out/gtdbtk.bac120.summary.tsv
        done < genome_manifest.tsv
    else
        mkdir -p gtdbtk_out
    fi
    : > gtdbtk_out/gtdbtk.ar53.summary.tsv

    cp gtdbtk_out/gtdbtk.bac120.summary.tsv gtdbtk.bac120.summary.tsv 2>/dev/null || : > gtdbtk.bac120.summary.tsv
    cp gtdbtk_out/gtdbtk.ar53.summary.tsv  gtdbtk.ar53.summary.tsv  2>/dev/null || : > gtdbtk.ar53.summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: 2.7.2
    END_VERSIONS
    """
}
