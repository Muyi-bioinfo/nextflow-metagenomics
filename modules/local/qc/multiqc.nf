// ============================================================================
// MULTIQC — 汇总 QC 报告
//
// 职责: 将各工具的 QC 输出汇总为单一 HTML 报告。
//
// Phase 3 输入源: FastQC (zip), fastp (json), Bowtie2 宿主去除 (log)
// 后续 Phase 追加: Kraken2, QUAST, mag_*.tsv 等 —— 通过向 ch_multiqc_files
// 追加通道实现, 追加新文件源时本模块无需改动。
//
// Phase 15: 增加 multiqc_config 输入 (multiqc_config.yaml, 经 -c 显式传入),
// 用 MultiQC custom content 表格式解析 mag_qc / mag_taxonomy / mag_metadata /
// mag_functional_annotation / mag_abundance 五张 TSV。另在暂存副本上清洗
// QUAST report.tsv 的 predicted genes 占位行 (见 script 内注释)。
//
// 注: 本 process 汇总的是整个 batch 的所有样本, 因此不带 meta —— 输入是
// collect() 后的文件集合。
// ============================================================================

process MULTIQC {
    label 'process_low'

    conda "bioconda::multiqc=1.35"
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'

    publishDir { "${params.outdir}/${params.batch_id}/99_multiqc" }, mode: 'copy'

    input:
    path multiqc_files, stageAs: "qc_inputs/*"
    path multiqc_config

    output:
    path "multiqc_report.html",       emit: report
    // 数据目录名由 --filename 派生 (multiqc_report.html -> multiqc_report_data),
    // 而非固定的 multiqc_data
    path "multiqc_report_data",       emit: data
    path "versions.yml",              emit: versions

    script:
    """
    # ---- Phase 15 上游数据清洗 (MultiQC 兼容性) ----
    # 1) Nextflow 默认以符号链接暂存目录输入 (如 QUAST 输出目录), 而 find
    #    不下钻链接目录 —— 先把暂存区里的链接目录替换为真实副本, 同时避免
    #    原地修改穿透链接污染上游 work 目录 (副本修改不影响原文件)。
    find qc_inputs -type l -print0 | while IFS= read -r -d '' l; do
        [ -d "\$l" ] || continue
        target=\$(readlink -f "\$l")
        rm "\$l"
        cp -r "\$target" "\$l"
    done

    # 2) QUAST 未做基因预测时 report.tsv 的 "# predicted genes (>= N bp)"
    #    行为 "-" —— MultiQC 的 quast 模块 (1.21 与 1.35 源码同段) 对这些
    #    值做字符串减法, 直接 TypeError 崩掉整个 quast 模块。删除"全部为 -"
    #    的占位行: 有真实数值 (做过基因预测) 的行保留, 不伪造数据; 只改
    #    qc_inputs/ 暂存副本, 发布目录 (05_assembly/quast/) 的原始文件不受影响。
    find qc_inputs -name "report.tsv" -print0 | while IFS= read -r -d '' f; do
        awk -F'\\t' 'BEGIN{OFS=FS} {keep=1; if (\$1 ~ /^# predicted genes/) {keep=0; for(i=2;i<=NF;i++) if(\$i != "-") {keep=1; break}} if(keep) print}' "\$f" > "\$f.tmp" && mv "\$f.tmp" "\$f"
    done

    multiqc \\
        --force \\
        --title "nextflow-metagenomics — ${params.batch_id}" \\
        --filename multiqc_report.html \\
        -c ${multiqc_config} \\
        qc_inputs/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed 's/multiqc, version //' )
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_report_data

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: 1.35
    END_VERSIONS
    """
}
