// ============================================================================
// MAG_QC_SUMMARY — 汇总每个 MAG 的 CheckM2 结果 (Phase 8)
// ============================================================================

process MAG_QC_SUMMARY {
    label 'process_single'

    conda "conda-forge::python=3.12"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/08_mag_qc" }, mode: 'copy'

    input:
    path qc_reports

    output:
    path "mag_qc.tsv", emit: tsv
    path "versions.yml", emit: versions

    script:
    def report_list = qc_reports instanceof List ? qc_reports.collect { "'${it}'" }.join(', ') : "'${qc_reports}'"
    """
    #!/usr/bin/env python3
    import csv
    from pathlib import Path

    reports = sorted([Path(p) for p in [${report_list}]])
    columns = ["meta_id", "mag_id", "Name", "Completeness", "Contamination"]
    with open("mag_qc.tsv", "w", newline="") as out:
        writer = csv.writer(out, delimiter="\\t")
        writer.writerow(columns)
        for report in reports:
            with report.open(newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\\t")
                for row in reader:
                    writer.writerow([row.get(c, "") for c in columns])

    with open("versions.yml", "w") as v:
        v.write('"${task.process}":\\n    python: "3.12"\\n')
    """

    stub:
    """
    printf "meta_id\\tmag_id\\tName\\tCompleteness\\tContamination\\n" > mag_qc.tsv
    for report in ${qc_reports}; do tail -n +2 "\$report" >> mag_qc.tsv; done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
