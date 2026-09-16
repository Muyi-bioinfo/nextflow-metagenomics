// ============================================================================
// DREP — MAG 去冗余 (Phase 9)
//
// dRep dereplicate 是**集合级批处理**工具: 本 process 把本阶段全部 qualified
// MAG 聚合为一个基因组目录, 一次调用完成 ANI 聚类 (MASH 预筛 + FastANI
// 精算, 两者随 dRep 环境自带, 无外部参考库依赖)。每簇按打分 (Phase 8
// CheckM2 的 completeness/contamination) 选一个代表基因组, 输出
// dereplicated_genomes/ (代表 FASTA) 与 data_tables/ (聚类原始数据表)。
//
// 输入:  genomeInfo.csv (genome,completeness,contamination —— 由 Phase 8 的
//        mag_qc.tsv 经 DREP 子工作流 join qualified MAG 后映射而来)
//        + 输入 MAG FASTA 路径列表 (toSortedList 聚合)。每个 MAG 以其 FASTA
//        文件名 (即 <mag_id>.fa, Phase 7 稳定 MAG ID 命名保证) 软链入
//        genomes/, dRep 输出表的 genome 列即 <mag_id>, 身份映射不丢失。
// 输出:  dereplicated_genomes/ + data_tables/ (dRep 原始输出, 移到 task
//        根目录按目录名发布, 保留不删)。
//
// 注意: 不用 val 通道传 tuple 列表 —— Nextflow 的 collect 聚合会把 tuple
// 扁平化为 ArrayBag, process 内无法按三元组还原 (实测发现)。
// ============================================================================

process DREP {
    tag "batch: ${params.batch_id}"
    label 'process_high'

    // Phase 17 审计: 对齐 environment.yml 的 3.7.1 (原 3.4.3 为 2021 年版本,
    // 与 env 的 3.7.1 不一致且无记录理由; tag 3.7.1--pyhdfd78af_0 已在 quay.io 验证存在)
    conda "bioconda::drep=3.7.1"
    container 'quay.io/biocontainers/drep:3.7.1--pyhdfd78af_0'

    // Phase 7 教训 —— publishDir 路径末尾与 pattern 都带目录名会产生
    // mags/mags/ 式双层嵌套, 这里 publishDir 路径不含 pattern 目录名。
    publishDir { "${params.outdir}/${params.batch_id}/09_dereplication" },
        mode: 'copy', pattern: "dereplicated_genomes"
    publishDir { "${params.outdir}/${params.batch_id}/09_dereplication" },
        mode: 'copy', pattern: "data_tables"
    publishDir { "${params.outdir}/${params.batch_id}/09_dereplication" },
        mode: 'copy', pattern: "drep_log"

    input:
    path(genome_info)
    path(fastas)

    output:
    path "dereplicated_genomes", type: 'dir', emit: catalog
    path "data_tables",          type: 'dir', emit: clusters
    path "versions.yml",          emit: versions

    script:
    def fasta_list = fastas instanceof List ? fastas : [fastas]
    def staging = fasta_list.collect { f -> "ln -s ${f} genomes/${f.getFileName()}" }.join('\n')
    """
    mkdir -p genomes

    # 同名 MAG ID 会在软链时互相覆盖, 且 dRep 以 genome 名为标识 —— 防御性
    # 检查 (Phase 7 稳定 MAG ID 保证唯一, 重名让成员表映射失效)
    dupes=\$(for f in ${fastas}; do basename "\$f"; done | sort | uniq -d)
    if [ -n "\$dupes" ]; then
        echo "ERROR: duplicate MAG FASTA basename — dRep genome names must be unique: \$dupes" >&2
        exit 1
    fi

    # 每个 MAG 以其 FASTA 文件名 (即 <mag_id>.fa) 软链入同一基因组目录
    ${staging}

    # 重试安全: dRep 不支持覆盖已有输出目录
    rm -rf drep_out

    dRep dereplicate drep_out \\
        -g genomes/*.fa \\
        --genomeInfo genomeInfo.csv \\
        -p ${task.cpus} \\
        ${params.drep_args}

    # 关键输出移到 task 根目录, 便于按目录名发布 (drep_out 其余内容如
    # figures/ 不发布, 数据表与代表 FASTA 保留)
    mv drep_out/dereplicated_genomes .
    mv drep_out/data_tables .
    mv drep_out/log drep_log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drep: 3.4.3
    END_VERSIONS
    """

    stub:
    """
    mkdir -p drep_out/dereplicated_genomes drep_out/data_tables drep_out/log

    # stub: 不真实调用 dRep —— 全部 MAG 聚为一簇, 字典序第一个 MAG 为簇代表
    # (非平凡成员映射: 覆盖成员表解析与 Phase 10 只分类代表集的接线)
    : > genomes.txt
    for f in ${fastas}; do basename "\$f" .fa >> genomes.txt; done
    sort -u genomes.txt -o genomes.txt

    if [ -s genomes.txt ]; then
        rep="\$(head -n 1 genomes.txt)"

        printf 'genome,secondary_cluster\\n' > drep_out/data_tables/Cdb.csv
        printf 'genome,primary_cluster,secondary_cluster\\n' > drep_out/data_tables/Wdb.csv
        while read -r g; do
            [ -n "\$g" ] || continue
            printf '%s,%s\\n' "\$g" 1 >> drep_out/data_tables/Cdb.csv
            printf '%s,%s,%s\\n' "\$g" 1 1 >> drep_out/data_tables/Wdb.csv
        done < genomes.txt

        for f in ${fastas}; do
            b="\$(basename "\$f" .fa)"
            if [ "\$b" = "\$rep" ]; then cp "\$f" drep_out/dereplicated_genomes/; fi
        done
    fi

    mv drep_out/dereplicated_genomes .
    mv drep_out/data_tables .
    mv drep_out/log drep_log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drep: 3.4.3
    END_VERSIONS
    """
}
