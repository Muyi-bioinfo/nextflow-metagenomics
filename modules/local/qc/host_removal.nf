// ============================================================================
// HOST_REMOVAL — Bowtie2 宿主序列去除
//
// 职责: 将清洁 reads 比对到宿主基因组, 保留未比对上的 (= 非宿主) reads。
//
// 输入:  tuple(meta, reads)   reads 为列表: SE=[R1], PE=[R1,R2]
//        path(index)          Bowtie2 索引文件集合
// 输出:  tuple(meta, reads)   SE=[nonhost_R1], PE=[nonhost_R1, nonhost_R2]
//        Bowtie2 比对日志 (MultiQC 可解析) + samtools flagstat/idxstats 统计
//
// 索引处理: 索引前缀在运行时由索引文件名推导 (见下方 script), 而非硬编码。
// 用户仅需提供 --host_index <前缀>, 例如 --host_index /data/ref/GRCh38
// (对应 GRCh38.1.bt2, GRCh38.2.bt2, ... GRCh38.rev.1.bt2 等文件)。
//
// 关于宿主 BAM: 默认不发布 —— 宿主 reads 对下游分析无用且体积巨大。
// 但保留它以计算 flagstat/idxstats, 这些统计是判断宿主污染程度的关键指标。
// ============================================================================

process HOST_REMOVAL {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::bowtie2=2.5.5 bioconda::samtools=1.24"
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'

    publishDir { "${params.outdir}/${params.batch_id}/02_host_removal/logs" },  mode: 'copy', pattern: "*.{log,flagstat,idxstats}"
    publishDir { "${params.outdir}/${params.batch_id}/02_host_removal/reads" }, mode: 'copy', pattern: "*_nonhost_R*.fastq.gz", enabled: params.save_host_removed

    input:
    tuple val(meta), path(reads)
    path index

    output:
    // glob 自动匹配 SE=1 个文件 / PE=2 个文件
    tuple val(meta), path("${meta.id}_nonhost_R*.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.bowtie2.log"),         emit: log
    tuple val(meta), path("${meta.id}.flagstat"),            emit: flagstat
    tuple val(meta), path("${meta.id}.idxstats"),            emit: idxstats
    path "versions.yml",                                     emit: versions

    script:
    // ------------------------------------------------------------------
    // 从暂存的索引文件推导 Bowtie2 索引前缀。
    // 支持小索引 (.bt2) 与大基因组索引 (.bt2l); 以 .rev.1.bt2* 为锚点,
    // 因为该后缀在一套索引中唯一, 不会误匹配。
    // ------------------------------------------------------------------
    if (meta.single_end) {
        """
        INDEX=\$( find -L . -name "*.rev.1.bt2" -o -name "*.rev.1.bt2l" | head -n 1 | sed 's/\\.rev\\.1\\.bt2l\\?\$//' )

        if [ -z "\$INDEX" ]; then
            echo "ERROR: 在暂存文件中未找到 Bowtie2 索引 (*.rev.1.bt2 / *.rev.1.bt2l)。" >&2
            echo "       请检查 --host_index 前缀是否正确, 例如 --host_index /data/ref/GRCh38" >&2
            exit 1
        fi

        bowtie2 \\
            -x "\$INDEX" \\
            -U ${reads[0]} \\
            --threads ${task.cpus} \\
            --un-gz ${meta.id}_nonhost_R1.fastq.gz \\
            2> ${meta.id}.bowtie2.log \\
            | samtools sort -@ ${task.cpus} -O bam -o ${meta.id}.host.bam -

        samtools index -@ ${task.cpus} ${meta.id}.host.bam
        samtools flagstat -@ ${task.cpus} ${meta.id}.host.bam > ${meta.id}.flagstat
        samtools idxstats ${meta.id}.host.bam > ${meta.id}.idxstats

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bowtie2: \$( bowtie2 --version 2>&1 | head -n 1 | sed 's/^.*bowtie2-align-s version //; s/ .*\$//' )
            samtools: \$( samtools --version | head -n 1 | sed 's/samtools //' )
        END_VERSIONS
        """
    } else {
        """
        INDEX=\$( find -L . -name "*.rev.1.bt2" -o -name "*.rev.1.bt2l" | head -n 1 | sed 's/\\.rev\\.1\\.bt2l\\?\$//' )

        if [ -z "\$INDEX" ]; then
            echo "ERROR: 在暂存文件中未找到 Bowtie2 索引 (*.rev.1.bt2 / *.rev.1.bt2l)。" >&2
            echo "       请检查 --host_index 前缀是否正确, 例如 --host_index /data/ref/GRCh38" >&2
            exit 1
        fi

        # --un-conc-gz: 写出双端均未比对上宿主的 read 对 (% 由 bowtie2 替换为 1/2)
        bowtie2 \\
            -x "\$INDEX" \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            --threads ${task.cpus} \\
            --un-conc-gz ${meta.id}_nonhost_R%.fastq.gz \\
            2> ${meta.id}.bowtie2.log \\
            | samtools sort -@ ${task.cpus} -O bam -o ${meta.id}.host.bam -

        samtools index -@ ${task.cpus} ${meta.id}.host.bam
        samtools flagstat -@ ${task.cpus} ${meta.id}.host.bam > ${meta.id}.flagstat
        samtools idxstats ${meta.id}.host.bam > ${meta.id}.idxstats

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bowtie2: \$( bowtie2 --version 2>&1 | head -n 1 | sed 's/^.*bowtie2-align-s version //; s/ .*\$//' )
            samtools: \$( samtools --version | head -n 1 | sed 's/samtools //' )
        END_VERSIONS
        """
    }

    stub:
    if (meta.single_end) {
        """
        echo | gzip > ${meta.id}_nonhost_R1.fastq.gz
        touch ${meta.id}.bowtie2.log ${meta.id}.flagstat ${meta.id}.idxstats

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bowtie2: 2.5.5
            samtools: 1.24
        END_VERSIONS
        """
    } else {
        """
        echo | gzip > ${meta.id}_nonhost_R1.fastq.gz
        echo | gzip > ${meta.id}_nonhost_R2.fastq.gz
        touch ${meta.id}.bowtie2.log ${meta.id}.flagstat ${meta.id}.idxstats

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bowtie2: 2.5.5
            samtools: 1.24
        END_VERSIONS
        """
    }
}
