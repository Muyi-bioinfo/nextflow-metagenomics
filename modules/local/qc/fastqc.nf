// ============================================================================
// FASTQC — 原始 reads 质量评估
//
// 职责: 仅评估质量, 不做任何过滤或修剪 (过滤由 FASTP 负责)。
//
// 输入:  tuple(meta, reads)   reads 为列表: SE=[R1], PE=[R1,R2]
// 输出:  HTML 报告 + ZIP 存档 (ZIP 供 MultiQC 汇总)
// ============================================================================

process FASTQC {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::fastqc=0.12.1"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    publishDir { "${params.outdir}/${params.batch_id}/01_qc/fastqc" }, mode: 'copy', pattern: "*.{html,zip}"

    input:
    // stageAs: 将输入暂存到子目录, 避免与下方重命名的目标文件名冲突 ——
    // 若输入本身已叫 <sample>_R1.fastq.gz, 同名 ln 会创建自指向的坏链接。
    tuple val(meta), path(reads, stageAs: 'input/*')

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip
    path "versions.yml",             emit: versions

    script:
    // 重命名为 <sample>_R{1[,2]}: FastQC 以输入文件名作为报告标题与 zip 内的
    // 样本名, 统一命名可确保 MultiQC 中的样本名与 meta.id 一致。
    if (meta.single_end) {
        """
        ln -s ${reads[0]} ${meta.id}_R1.fastq.gz

        fastqc \\
            --threads ${task.cpus} \\
            --quiet \\
            ${meta.id}_R1.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastqc: \$( fastqc --version | sed 's/FastQC v//' )
        END_VERSIONS
        """
    } else {
        """
        ln -s ${reads[0]} ${meta.id}_R1.fastq.gz
        ln -s ${reads[1]} ${meta.id}_R2.fastq.gz

        fastqc \\
            --threads ${task.cpus} \\
            --quiet \\
            ${meta.id}_R1.fastq.gz \\
            ${meta.id}_R2.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastqc: \$( fastqc --version | sed 's/FastQC v//' )
        END_VERSIONS
        """
    }

    stub:
    if (meta.single_end) {
        """
        touch ${meta.id}_R1_fastqc.html
        touch ${meta.id}_R1_fastqc.zip

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastqc: 0.12.1
        END_VERSIONS
        """
    } else {
        """
        touch ${meta.id}_R1_fastqc.html ${meta.id}_R2_fastqc.html
        touch ${meta.id}_R1_fastqc.zip  ${meta.id}_R2_fastqc.zip

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastqc: 0.12.1
        END_VERSIONS
        """
    }
}
