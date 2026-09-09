// ============================================================================
// FASTP — 接头修剪与质量过滤
//
// 职责: 生成清洁 reads, 并输出结构化 QC 指标。
//
// 输入:  tuple(meta, reads)   reads 为列表: SE=[R1], PE=[R1,R2]
// 输出:  tuple(meta, reads)   SE=[clean_R1], PE=[clean_R1, clean_R2]
//        JSON (机器可读, 供 MultiQC 与指标汇总) + HTML (人工查看)
//
// JSON 中记录的指标: raw reads / filtered reads / Q20 / Q30 / adapter /
// duplication / GC — 由 bin/parse_fastp_json.py 汇总为 TSV。
//
// 关键: 本 process 是全流程唯一的 clean reads 来源。read-based 分支与
// assembly 分支共用同一份输出, 不得各自重新生成 (见 subworkflows/local/preprocessing.nf)。
// ============================================================================

process FASTP {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::fastp=1.3.6"
    container 'quay.io/biocontainers/fastp:1.3.6--h43da1c4_0'

    publishDir { "${params.outdir}/${params.batch_id}/01_qc/fastp" }, mode: 'copy', pattern: "*.{json,html}"
    publishDir { "${params.outdir}/${params.batch_id}/01_qc/fastp/clean_reads" }, mode: 'copy', pattern: "*.fastp.fastq.gz", enabled: params.save_trimmed

    input:
    tuple val(meta), path(reads)

    output:
    // glob 自动匹配 SE=1 个文件 / PE=2 个文件, 下游统一用 reads[0] / reads[1] 访问
    tuple val(meta), path("*.fastp.fastq.gz"),       emit: reads
    tuple val(meta), path("${meta.id}.fastp.json"),  emit: json
    tuple val(meta), path("${meta.id}.fastp.html"),  emit: html
    tuple val(meta), path("${meta.id}.fastp.log"),   emit: log
    path "versions.yml",                             emit: versions

    script:
    def dedup_arg = params.fastp_dedup ? '--dedup' : ''
    if (meta.single_end) {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --out1 ${meta.id}_R1.fastp.fastq.gz \\
            --json ${meta.id}.fastp.json \\
            --html ${meta.id}.fastp.html \\
            --report_title "${meta.id} — fastp report" \\
            --qualified_quality_phred ${params.fastp_qualified_quality} \\
            --unqualified_percent_limit ${params.fastp_unqualified_percent} \\
            --length_required ${params.fastp_min_length} \\
            --cut_front \\
            --cut_tail \\
            --cut_mean_quality ${params.fastp_cut_mean_quality} \\
            ${dedup_arg} \\
            --thread ${task.cpus} \\
            2> >( tee ${meta.id}.fastp.log >&2 )

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$( fastp --version 2>&1 | sed -e 's/fastp //g' )
        END_VERSIONS
        """
    } else {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --in2 ${reads[1]} \\
            --out1 ${meta.id}_R1.fastp.fastq.gz \\
            --out2 ${meta.id}_R2.fastp.fastq.gz \\
            --json ${meta.id}.fastp.json \\
            --html ${meta.id}.fastp.html \\
            --report_title "${meta.id} — fastp report" \\
            --detect_adapter_for_pe \\
            --qualified_quality_phred ${params.fastp_qualified_quality} \\
            --unqualified_percent_limit ${params.fastp_unqualified_percent} \\
            --length_required ${params.fastp_min_length} \\
            --cut_front \\
            --cut_tail \\
            --cut_mean_quality ${params.fastp_cut_mean_quality} \\
            ${dedup_arg} \\
            --thread ${task.cpus} \\
            2> >( tee ${meta.id}.fastp.log >&2 )

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$( fastp --version 2>&1 | sed -e 's/fastp //g' )
        END_VERSIONS
        """
    }

    stub:
    if (meta.single_end) {
        """
        echo | gzip > ${meta.id}_R1.fastp.fastq.gz
        touch ${meta.id}.fastp.json ${meta.id}.fastp.html ${meta.id}.fastp.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: 1.3.6
        END_VERSIONS
        """
    } else {
        """
        echo | gzip > ${meta.id}_R1.fastp.fastq.gz
        echo | gzip > ${meta.id}_R2.fastp.fastq.gz
        touch ${meta.id}.fastp.json ${meta.id}.fastp.html ${meta.id}.fastp.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: 1.3.6
        END_VERSIONS
        """
    }
}
