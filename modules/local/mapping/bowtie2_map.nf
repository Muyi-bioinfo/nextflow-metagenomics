// ============================================================================
// BOWTIE2_MAP — 将样本 clean reads 比对回组装 contigs
//
// 职责: 只做比对, 产出未排序 BAM。排序/索引/统计交给 samtools_sort.nf,
//       覆盖度计算交给 coverage.nf —— 三者职责不重叠。
//
// 输入:  tuple(meta, reads, index)
//          reads  SE=[R1] / PE=[R1,R2]   来自 PREPROCESSING (clean reads)
//          index  BOWTIE2_BUILD 产出的索引目录
//          meta   由 subworkflows/local/mapping.nf 配对生成, 同时携带
//                 「哪个样本」(sample / single_end) 与「比对到哪套组装」
//                 (assembly_id / assembler / assembly_mode / samples)
// 输出:  tuple(meta, bam)   未排序 BAM
//        bowtie2 日志 (MultiQC 可解析, 含总体比对率)
//
// ─── 与 HOST_REMOVAL 的区别 ───────────────────────────────────────────────
// 同样是 bowtie2, 但目的相反:
//   HOST_REMOVAL  比对到**宿主**, 保留**未**比对上的 reads (丢弃比对结果)
//   BOWTIE2_MAP   比对到**自己的 contigs**, 保留**比对结果本身** (BAM)
//                 —— BAM 里的深度信息正是分箱所需的覆盖度来源
//
// ─── 为什么保留未比对 reads ───────────────────────────────────────────────
// 默认不加 --no-unal: 保留未比对记录后 samtools flagstat 才能给出真实的
// 比对率分母, 这是判断「组装是否代表了样本主体」的关键指标。真实数据若在意
// BAM 体积, 可用 --bowtie2_map_args '--no-unal' 显式丢弃 (覆盖度计算不受
// 影响 —— jgi_summarize_bam_contig_depths 本就忽略未比对记录)。
// ============================================================================

process BOWTIE2_MAP {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::bowtie2=2.5.5 bioconda::samtools=1.24"
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'

    publishDir { "${params.outdir}/${params.batch_id}/06_mapping/logs" },
        mode: 'copy', pattern: "*.bowtie2.log"

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path("${meta.id}.bam"),         emit: bam
    tuple val(meta), path("${meta.id}.bowtie2.log"), emit: log
    path "versions.yml",                             emit: versions

    script:
    // 索引前缀在运行时从暂存文件反推, 不做约定式拼接 (同 host_removal.nf)
    def input_arg = meta.single_end ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    set -o pipefail

    INDEX=\$( find -L . -name "*.rev.1.bt2" -o -name "*.rev.1.bt2l" | head -n 1 | sed 's/\\.rev\\.1\\.bt2l\\?\$//' )

    if [ -z "\$INDEX" ]; then
        echo "ERROR: 在暂存文件中未找到 Bowtie2 索引 (*.rev.1.bt2 / *.rev.1.bt2l)。" >&2
        echo "       该索引由 BOWTIE2_BUILD 产出, 请检查上游任务是否成功。" >&2
        exit 1
    fi

    bowtie2 \\
        -x "\$INDEX" \\
        ${input_arg} \\
        --threads ${task.cpus} \\
        ${params.bowtie2_map_args} \\
        2> ${meta.id}.bowtie2.log \\
        | samtools view -@ ${task.cpus} -b -o ${meta.id}.bam -

    if [ ! -s ${meta.id}.bam ]; then
        echo "ERROR: 未产出 BAM (${meta.id}.bam 为空)。bowtie2 日志:" >&2
        cat ${meta.id}.bowtie2.log >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$( bowtie2 --version 2>&1 | head -n 1 | sed 's/^.*bowtie2-align-s version //; s/ .*\$//' )
        samtools: \$( samtools --version | head -n 1 | sed 's/samtools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.bam ${meta.id}.bowtie2.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: 2.5.5
        samtools: 1.24
    END_VERSIONS
    """
}
