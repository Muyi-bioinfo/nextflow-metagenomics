// ============================================================================
// DREP_SUMMARY — dRep 聚类输出 → 成员表 (Phase 9)
//
// 解析 dRep 的 Cdb.csv (簇归属) 与 dereplicated_genomes/ (簇代表), 结合
// manifest (mag_id → sample) 生成覆盖**全部 qualified MAG** 的成员表:
//     mag_id \t sample \t representative_mag_id   (无表头)
// 代表 MAG 的行 mag_id == representative_mag_id。映射不一致 (重复/缺失/
// 一簇多代表) 由 bin/derep_summary.py 直接报错, 不静默回退。
//
// 成员表是 Phase 10 分类回填与 Phase 9 代表集恢复 (join 回 meta/FASTA)
// 的共同依据。
// ============================================================================

process DREP_SUMMARY {
    tag "batch: ${params.batch_id}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/09_dereplication" },
        mode: 'copy', pattern: "mag_membership.tsv"

    input:
    path manifest        // genome_manifest.tsv: mag_id \t sample
    path data_tables     // dRep data_tables/ 目录
    path dereplicated    // dRep dereplicated_genomes/ 目录
    path parser          // bin/derep_summary.py 作为显式输入暂存, 使本 process 可从任意入口复用

    output:
    path "mag_membership.tsv", emit: membership
    path "versions.yml",       emit: versions

    script:
    """
    python3 ${parser} \\
        --manifest ${manifest} \\
        --cdb ${data_tables}/Cdb.csv \\
        --dereplicated ${dereplicated} \\
        --output mag_membership.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """

    stub:
    """
    # stub 与真实脚本一致: 输入 (manifest/Cdb/代表目录) 在 stub 模式下同样由
    # DREP 的 stub 产出, 直接跑真实解析逻辑, 成员表路径得到完整覆盖
    python3 ${parser} \\
        --manifest ${manifest} \\
        --cdb ${data_tables}/Cdb.csv \\
        --dereplicated ${dereplicated} \\
        --output mag_membership.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
