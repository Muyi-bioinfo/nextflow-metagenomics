// ============================================================================
// HUMANN — 功能谱分析 (gene family / pathway abundance / pathway coverage)
//
// ─── 与 Kraken2/Bracken 的关系: 并行, 不是下游 ────────────────────────────
// HUMAnN 直接消费 clean reads, 与 Kraken2 分支同级并行。它内部自带 prescreen
// (MetaPhlAn) 来决定检索哪些物种的基因组, 不需要 Kraken2 的输出。
// 正确:  clean reads ─┬─→ Kraken2 → Bracken
//                     └─→ HUMAnN
// 错误:  Kraken2 → Bracken → HUMAnN
//
// (可选) --taxonomic_profile 能跳过 prescreen, 但那需要 MetaPhlAn 格式的 profile,
// 而非 Kraken2 报告 —— 二者格式不兼容, 不能直接串联。
//
// ─── 双端输入的处理 ───────────────────────────────────────────────────────
// HUMAnN 只接受**单个**输入文件, 没有 -1/-2 参数。双端数据的官方做法是把 R1 与
// R2 拼接成一个文件 —— HUMAnN 按 read 独立比对, 不使用配对信息。
// 这里显式 cat 两个 gz (gzip 格式可直接拼接), 而非静默只用 R1 —— 后者会悄悄
// 丢掉一半数据, 且从输出上看不出来。
//
// ─── 三个输出的区别 ───────────────────────────────────────────────────────
//   genefamilies    UniRef 基因家族丰度 (RPK) —— 最细粒度
//   pathabundance   MetaCyc 通路丰度       —— 由基因家族聚合而来
//   pathcoverage    通路覆盖度 (0-1)       —— 通路的"完整程度", 与丰度是两个维度:
//                                             低丰度但高覆盖 = 存在但表达少
//                                             高丰度但低覆盖 = 可能是少数基因的假阳性
//
// ─── 数据库参数化 ─────────────────────────────────────────────────────────
// HUMAnN 需要两个独立数据库:
//   nucleotide (ChocoPhlAn) — params.humann_nucleotide_db
//   protein    (UniRef)     — params.humann_protein_db
// 另有 MetaPhlAn 数据库 (params.metaphlan_db) 供 prescreen 使用。
// 三者缺一不可; 未提供时本 process 不会被调用 (见 read_based.nf 的条件判断)。
//
// ─── 版本: 3.9 —— 与数据库世代绑定, 不可随意降级 ──────────────────────────
// 本模块曾声明 humann 3.0.1, 与 environment.yml 的 3.9 不一致。二者不是"新旧
// 之别"而是**数据库不兼容**: 3.0.x 走 MetaPhlAn 3 (mpa_v30 标记库 +
// ChocoPhlAn v296), 3.6+ 走 MetaPhlAn 4 (CHOCOPhlAnSGB 标记库 + ChocoPhlAn v31)。
// 同一个 --metaphlan_db / --humann_nucleotide_db 路径喂不了两代。
// 若模块与 environment.yml 各行其是, `-profile conda` 与 `-profile standard`
// 会要求两套完全不同的数据库, 且报错发生在下载 40 GB 之后。故统一到 3.9。
// (3.9 亦是 python 3.12 下唯一可用的 build, 见 environment.yml 顶部说明。)
// ============================================================================

process HUMANN {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::humann=3.9"
    container 'quay.io/biocontainers/humann:3.9--py312hdfd78af_0'

    publishDir { "${params.outdir}/${params.batch_id}/04_function/humann" }, mode: 'copy',
        pattern: "*.{tsv,log}"

    input:
    tuple val(meta), path(reads)
    path nucleotide_db
    path protein_db
    path metaphlan_db

    output:
    tuple val(meta), path("${meta.id}_genefamilies.tsv"),  emit: genefamilies
    tuple val(meta), path("${meta.id}_pathabundance.tsv"), emit: pathabundance
    tuple val(meta), path("${meta.id}_pathcoverage.tsv"),  emit: pathcoverage
    tuple val(meta), path("${meta.id}.humann.log"),        emit: log
    path "versions.yml",                                   emit: versions

    script:
    // MetaPhlAn 数据库索引名: 目录中 *.pkl 的基名。
    // MetaPhlAn 4 (humann 3.9 所需) 的索引名形如 mpa_v<版本>_CHOCOPhlAnSGB_<日期>;
    // MetaPhlAn 3 则是 mpa_v30_CHOCOPhlAn_201901 —— 后者与本模块的 3.9 不兼容。
    def metaphlan_opts = "--bowtie2db ${metaphlan_db}"
    """
    # ------------------------------------------------------------------
    # 双端: 拼接为单文件。HUMAnN 无 -1/-2, 按 read 独立比对。
    # gzip 成员可直接串联, 解压后等价于两个文件依次拼接。
    # ------------------------------------------------------------------
    ${ meta.single_end
        ? "cat ${reads[0]} > ${meta.id}.humann_input.fastq.gz"
        : "cat ${reads[0]} ${reads[1]} > ${meta.id}.humann_input.fastq.gz" }

    # MetaPhlAn 索引名从 .pkl 推导, 不硬编码版本 —— 数据库更新后无需改模块
    MPA_INDEX=\$( find -L ${metaphlan_db} -name "*.pkl" | head -n 1 | xargs -r basename | sed 's/\\.pkl\$//' )

    if [ -z "\$MPA_INDEX" ]; then
        echo "ERROR: 在 ${metaphlan_db} 中未找到 MetaPhlAn 索引 (*.pkl)。" >&2
        echo "       --metaphlan_db 应指向 metaphlan --install 生成的数据库目录。" >&2
        exit 1
    fi

    humann \\
        --input ${meta.id}.humann_input.fastq.gz \\
        --output . \\
        --output-basename ${meta.id} \\
        --threads ${task.cpus} \\
        --nucleotide-database ${nucleotide_db} \\
        --protein-database ${protein_db} \\
        --metaphlan-options "${metaphlan_opts} --index \$MPA_INDEX" \\
        --o-log ${meta.id}.humann.log \\
        --remove-temp-output \\
        ${params.humann_args}

    # 拼接产物体积等同输入, 不保留
    rm -f ${meta.id}.humann_input.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: \$( humann --version 2>&1 | sed 's/^humann v//' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_genefamilies.tsv \\
          ${meta.id}_pathabundance.tsv \\
          ${meta.id}_pathcoverage.tsv \\
          ${meta.id}.humann.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: 3.9
    END_VERSIONS
    """
}
