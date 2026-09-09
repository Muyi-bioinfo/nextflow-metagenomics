// ============================================================================
// CHECKM2 — MAG 质量控制 (Phase 8)
//
// 输入:  tuple(meta, mag_id, mag_fasta)
// 输出:  tuple(meta, mag_id, mag_fasta, mag_qc.tsv)
// ============================================================================

process CHECKM2 {
    tag "${meta.id}.${mag_id}"
    label 'process_high'

    // Phase 17 审计: 对齐 environment.yml 的 1.1.0 (build _1, python>3.12);
    // 原 1.0.2--pyhdfd78af_0 tag 在 quay.io 上不存在 (真 1.0.2 build 是
    // pyh7cba7a3_0), 容器模式拉取必然失败。
    conda "bioconda::checkm2=1.1.0"
    container 'quay.io/biocontainers/checkm2:1.1.0--pyh7e72e81_1'

    publishDir { "${params.outdir}/${params.batch_id}/08_mag_qc/checkm2" },
        mode: 'copy', pattern: "*.qc.tsv"

    input:
    tuple val(meta), val(mag_id), path(mag_fasta)

    output:
    tuple val(meta), val(mag_id), path("${mag_id}.fa"), path("${mag_id}.qc.tsv"), emit: results
    path "versions.yml", emit: versions

    script:
    def db_arg = "--database ${params.checkm2_db}"
    """
    # Use the stable MAG ID as the CheckM2 input basename so its Name column
    # remains directly mappable to the upstream MAG record.
    mkdir -p checkm2_input checkm2_out
    cp ${mag_fasta} checkm2_input/${mag_id}.fa

    checkm2 predict \\
        --input checkm2_input \\
        --output-directory checkm2_out \\
        --threads ${task.cpus} \\
        ${db_arg} \\
        ${params.checkm2_args}

    python3 - <<'PY'
    import csv
    from pathlib import Path

    reports = sorted(Path("checkm2_out").glob("quality_report.tsv"))
    if not reports:
        reports = sorted(Path("checkm2_out").glob("*.tsv"))
    if not reports:
        raise SystemExit("CheckM2 did not produce a quality report TSV")

    with reports[0].open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\\t"))
    if not rows:
        raise SystemExit("CheckM2 quality report contains no MAG records")

    row = rows[0]
    def find_column(names):
        for name in names:
            for key in row:
                if key.strip().lower() == name:
                    return key
        return None

    name_col = find_column(["name"])
    comp_col = find_column(["completeness"])
    cont_col = find_column(["contamination"])
    if not comp_col or not cont_col:
        raise SystemExit("CheckM2 report lacks Completeness/Contamination columns")

    with open("${mag_id}.qc.tsv", "w", newline="") as out:
        writer = csv.writer(out, delimiter="\\t")
        writer.writerow(["meta_id", "mag_id", "Name", "Completeness", "Contamination"])
        writer.writerow([
            "${meta.id}", "${mag_id}", row.get(name_col, "${mag_id}"),
            row[comp_col], row[cont_col]
        ])
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: 1.0.2
    END_VERSIONS
    """

    stub:
    """
    # 上游 stub 产出的 FASTA 可能已以 ${mag_id}.fa 命名 (与目标同名), cp 同文件
    # 会报 "are the same file"; 已同名时跳过。仅影响 stub, 真实脚本不受影响。
    [ -e ${mag_id}.fa ] || cp ${mag_fasta} ${mag_id}.fa
    printf "meta_id\\tmag_id\\tName\\tCompleteness\\tContamination\\n" > ${mag_id}.qc.tsv
    printf "${meta.id}\\t${mag_id}\\t${mag_id}\\t75.0\\t5.0\\n" >> ${mag_id}.qc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: 1.0.2
    END_VERSIONS
    """
}
