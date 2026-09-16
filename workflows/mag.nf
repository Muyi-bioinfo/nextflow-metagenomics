// ============================================================================
// nextflow-metagenomics — 主工作流 (workflows/mag.nf)
//
// Phase 3: 预处理    Phase 4: read-based 分支    Phase 5: 组装    Phase 6: 比对    Phase 7: 分箱
//
// 已实现:
//   Phase 2  样本表解析与验证 → tuple(meta, reads)
//            meta.single_end = true/false; reads=[R1] (SE) 或 [R1,R2] (PE)
//   Phase 3  FastQC → fastp → Bowtie2 宿主去除 → clean reads → MultiQC
//   Phase 4  clean reads ─┬─→ Kraken2 → Bracken   (03_taxonomy)
//                         └─→ HUMAnN              (04_function)
//   Phase 5  clean reads ──→ MEGAHIT ∥ metaSPAdes → QUAST  (05_assembly)
//   Phase 6  clean reads + contigs → Bowtie2 → samtools → 覆盖度  (06_mapping)
//   Phase 7  contigs + depth → MetaBAT2 → split bins → MAG FASTA  (07_binning)
//   Phase 8  MAG FASTA → CheckM2 QC → qualified MAGs              (08_mag_qc)
//   Phase 9  qualified MAGs → dRep 去冗余 → 代表 MAG + 成员表     (09_dereplication)
//   Phase 10 代表 MAGs → GTDB-Tk → mag_taxonomy.tsv               (10_mag_taxonomy)
//   Phase 11 代表 MAGs → Prodigal → genes.fna / proteins.faa / GFF (11_gene_prediction)
//   Phase 12 代表 MAGs → DIAMOND ∥ eggNOG-mapper ∥ RGI → 键控汇总表   (12_annotation)
//   Phase 13 代表 MAGs + 排序 BAM → CoverM → 跨样本丰度矩阵           (13_abundance)
//   Phase 14  各 Phase 汇总表 → 整合结果集 (两张核心表 + 引用拷贝)     (14_integrated)
//
// read-based 分支 (Phase 4) 与 assembly 分支 (Phase 5) 是**并行**的:
// 二者都从同一个 PREPROCESSING.out.reads 通道取数据, 彼此没有依赖, 也不共享
// 中间产物。Nextflow 会同时调度两条分支, 不会重复执行上游任务。
//
// Phase 6 则**串在组装分支之后**: 它把 clean reads 比回 Phase 5 的 contigs,
// 产出 MetaBAT2 所需的覆盖度矩阵。
//
// Phase 7 则**串在 Phase 5 + 6 之后**: 它以 contigs (Phase 5) + depth (Phase 6)
// 为输入, MetaBAT2 分箱后产出带稳定 MAG ID 的独立 FASTA。
//
// Phase 8 (CheckM2 质控) 串在 Phase 7 之后; Phase 9 (dRep 去冗余)、Phase 10
// (GTDB-Tk 分类) 与 Phase 11 (Prodigal 基因预测) 串在其后 —— 顺序按最佳
// 实践为**先去冗余后分类后注释** (GTDB-Tk 与 Prodigal 都只跑代表 MAG,
// 冗余成员的分类由 dRep 成员表回填)。四者都消费每个 MAG 一条记录
// (tuple(meta, mag_id, mag_fasta)) 的通道。
// ============================================================================

include { PREPROCESSING } from '../subworkflows/local/preprocessing.nf'
include { READ_BASED    } from '../subworkflows/local/read_based.nf'
include { READ_BASED_MERGE } from '../subworkflows/local/read_based_merge.nf'
include { ASSEMBLY      } from '../subworkflows/local/assembly.nf'
include { MAPPING       } from '../subworkflows/local/mapping.nf'
include { BINNING       } from '../subworkflows/local/binning.nf'
include { MAG_QC        } from '../subworkflows/local/mag_qc.nf'
include { DEREPLICATION as DREP } from '../subworkflows/local/dereplication.nf'
include { TAXONOMY      } from '../subworkflows/local/taxonomy.nf'
include { GENE_PREDICTION } from '../subworkflows/local/gene_prediction.nf'
include { ANNOTATION      } from '../subworkflows/local/annotation.nf'
include { ABUNDANCE       } from '../subworkflows/local/abundance.nf'
include { INTEGRATION     } from '../subworkflows/local/integration.nf'
include { PLOTTING        } from '../subworkflows/local/plotting.nf'
include { MULTIQC       } from '../modules/local/qc/multiqc.nf'

/*
 * 样本表解析和验证
 */
process CHECK_SAMPLESHEET {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container 'quay.io/biocontainers/python:3.12'

    publishDir { "${params.outdir}/${params.batch_id}/00_metadata" }, mode: 'copy'

    input:
    path samplesheet
    // 验证脚本作为显式输入暂存 (而非依赖 PATH), 使本 process 可从任意入口复用
    path validator
    // launchDir 以字符串传入: Path 对象参与 task hash 时不稳定, 会导致
    // -resume 无法命中缓存
    val base_dir

    output:
    path "validated_samplesheet.csv", emit: csv

    script:
    """
    python3 ${validator} $samplesheet --base-dir $base_dir > validated_samplesheet.csv
    """
}

workflow mag {

    /*
     * Batch ID 与输出目录:
     * params.batch_id 在 nextflow.config 中解析。各 process 的 publishDir 写成
     * 闭包形式 { "${params.outdir}/${params.batch_id}/..." } —— 闭包延迟到 task
     * 提交时求值, 使 profile 与命令行对 outdir 的覆盖能够生效。
     *
     * batch directory 与 Nextflow work/ cache 职责不同:
     *   <outdir>/<batch_id>/  -> 业务结果
     *   work/                 -> Nextflow task cache (-resume 依赖此目录)
     */

    if (!params.input) {
        error "ERROR: --input parameter is required. Please provide a samplesheet CSV file."
    }

    log.info """
    ======================================================
     nextflow-metagenomics
     input:         ${params.input}
     batch:         ${params.batch_id}
     outdir:        ${params.outdir}/${params.batch_id}
     host_index:    ${params.host_index ?: '(未提供 — 将跳过宿主去除)'}
     kraken2_db:    ${params.kraken2_db ?: '(未提供 — 将跳过 Kraken2/Bracken)'}
     humann_db:     ${params.humann_db ?: params.humann_nucleotide_db ?: '(未提供 — 将跳过 HUMAnN)'}
     assembler:     ${params.skip_assembly ? '(已跳过 — --skip_assembly)' : params.assembler}
     assembly_mode: ${params.assembly_mode}
     threads:       ${params.threads}
     execution dir: ${workflow.launchDir}
    ======================================================
    """.stripIndent()

    // ---------------------------------------------------------------------
    // Phase 2: 输入验证与元数据解析
    // ---------------------------------------------------------------------
    ch_input_csv = Channel.fromPath(params.input, checkIfExists: true)
    ch_validator = file("${projectDir}/bin/check_samplesheet.py", checkIfExists: true)
    CHECK_SAMPLESHEET(ch_input_csv, ch_validator, workflow.launchDir.toString())

    ch_raw_reads = CHECK_SAMPLESHEET.out.csv
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [
                id:         row.sample,
                group:      row.group,
                batch:      row.batch,
                host:       row.host,
                single_end: row.single_end == 'True'
            ]
            // 单端: reads = [R1]; 双端: reads = [R1, R2]
            def reads = meta.single_end
                ? [ file(row.fastq_1, checkIfExists: true) ]
                : [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ]

            return tuple(meta, reads)
        }

    // ---------------------------------------------------------------------
    // Phase 3: 预处理
    // ---------------------------------------------------------------------
    PREPROCESSING(ch_raw_reads)

    // clean reads —— read-based 与 assembly 分支从此通道接入, 两者共用同一份
    // 文件 (通道可被多个下游 process 消费, 上游不重复执行)
    ch_clean_reads   = PREPROCESSING.out.reads
    ch_multiqc_files = PREPROCESSING.out.multiqc_files

    // ---------------------------------------------------------------------
    // Phase 4: read-based 分支 (Kraken2 → Bracken ∥ HUMAnN)
    //
    // 该子工作流内部是两条并行分支; 各工具缺数据库时自行 warn 并跳过, 不会
    // 影响其他分支, 也不会伪造产物。
    // ---------------------------------------------------------------------
    ch_bracken_merged       = Channel.empty()
    ch_beta_diversity       = Channel.empty()
    ch_pathabundance_merged = Channel.empty()

    if (!params.skip_read_based) {
        READ_BASED(ch_clean_reads)
        ch_multiqc_files = ch_multiqc_files.mix(READ_BASED.out.multiqc_files)

        // Phase 20: read-based 跨样本整合 (只消费逐样本 emit, 产出样本×taxa /
        // 样本×pathway 宽表矩阵 + β 多样性距离矩阵)。skip_kraken2 / skip_bracken /
        // skip_humann 各自跳过时, 对应逐样本通道为空 → 对应 merge process 不
        // 调度 (告警已由 read_based.nf 发出)。
        READ_BASED_MERGE(
            READ_BASED.out.bracken_abundance,
            READ_BASED.out.pathabundance
        )
        ch_bracken_merged       = READ_BASED_MERGE.out.bracken_merged
        ch_beta_diversity       = READ_BASED_MERGE.out.beta_diversity
        ch_pathabundance_merged = READ_BASED_MERGE.out.pathabundance_merged
    }
    else {
        log.warn "整个 read-based 分支已跳过 (--skip_read_based), 不产出物种分类与功能谱结果。"
    }

    // ---------------------------------------------------------------------
    // Phase 5: 组装分支 (MEGAHIT ∥ metaSPAdes → QUAST)
    //
    // 与 Phase 4 并行, 同样从 ch_clean_reads 接入。子工作流内部按
    // params.assembler 选择组装器 —— megahit / metaspades 是可选替代,
    // "both" 时并行各跑一次以作比较, 而非串联。
    // ---------------------------------------------------------------------
    ch_contigs = Channel.empty()
    ch_assembly_summary = Channel.empty()

    if (!params.skip_assembly) {
        ASSEMBLY(ch_clean_reads)

        // Phase 6 (reads 比对) 与 Phase 7 (MetaBAT2 分箱) 从此通道接入
        ch_contigs       = ASSEMBLY.out.contigs
        // 组装汇总表 —— Phase 14 只做 14_integrated/ 下的引用拷贝, 不重做
        // (--skip_quast 时该通道为空, 拷贝自然不调度)
        ch_assembly_summary = ASSEMBLY.out.assembly_summary
        ch_multiqc_files = ch_multiqc_files.mix(ASSEMBLY.out.multiqc_files)
    }
    else {
        log.warn "整个组装分支已跳过 (--skip_assembly), 不产出 contigs, 后续 MAG 相关分析也无从进行。"
    }

    // ---------------------------------------------------------------------
    // Phase 6: reads 比对与覆盖度 (Bowtie2 → samtools → contig depth)
    //
    // 与 Phase 4/5 不同, 本阶段**依赖**组装分支: 比对的参考序列就是 Phase 5
    // 产出的 contigs。--skip_assembly 时 ch_contigs 为空通道, 这里不会有任何
    // 任务被调度 (上面已就此告警, 不重复)。
    //
    // clean reads 在此第三次被消费 (read-based / 组装 / 比对), 三者共用同一份
    // 文件, 上游不重复执行。
    // ---------------------------------------------------------------------
    ch_depth = Channel.empty()

    if (!params.skip_mapping) {
        MAPPING(ch_contigs, ch_clean_reads)

        // Phase 7 (MetaBAT2 分箱) 从此通道接入 —— 与 ch_contigs 按
        // (组装单元 id, assembler) 配对即可
        ch_depth         = MAPPING.out.depth
        ch_multiqc_files = ch_multiqc_files.mix(MAPPING.out.multiqc_files)
    }
    else {
        log.warn "reads 比对已跳过 (--skip_mapping), 不产出 BAM 与覆盖度矩阵, 后续 MAG 分箱 (Phase 7) 也无从进行。"
    }

    // ---------------------------------------------------------------------
    // Phase 7: MAG 分箱 (MetaBAT2 → split bins → MAG FASTA)
    //
    // 与 Phase 6 串联: 需要 contigs (Phase 5) + depth (Phase 6)。
    // --skip_assembly 或 --skip_mapping 时两个上游通道至少有一个为空, 这里不会
    // 有任何任务被调度 (上面已就此告警, 不重复)。
    // ---------------------------------------------------------------------
    ch_mags        = Channel.empty()
    ch_bin_summary = Channel.empty()

    if (!params.skip_binning) {
        BINNING(ch_contigs, ch_depth)

        // Phase 8 (CheckM2 质控) 从此通道接入 —— 每个 MAG 一条记录
        ch_mags          = BINNING.out.mags
        ch_bin_summary   = BINNING.out.summary   // Phase 21 漏斗 raw bins
        ch_multiqc_files = ch_multiqc_files.mix(BINNING.out.multiqc_files)

        // Phase 8: CheckM2 MAG QC; qualified_mags retains BINNING.out.mags shape.
        MAG_QC(ch_mags)
        ch_qualified_mags = MAG_QC.out.qualified_mags
        ch_qc_table       = MAG_QC.out.qc_table
        ch_multiqc_files  = ch_multiqc_files.mix(MAG_QC.out.multiqc_files)  // Phase 15: mag_qc.tsv
    }
    else {
        log.warn "MAG 分箱已跳过 (--skip_binning), 不产出 MAG FASTA, 后续 CheckM2 / GTDB-Tk / dRep 等 MAG 相关分析也无从进行。"
        ch_qualified_mags = Channel.empty()
        ch_qc_table       = Channel.empty()
    }

    // ---------------------------------------------------------------------
    // Phase 9: MAG 去冗余 (dRep dereplicate → 代表 MAG + 成员表)
    //
    // 集合级批处理: 全部 qualified MAG 聚合为一个基因组目录, dRep 一次调用
    // 完成 ANI 聚类, 每簇选一个代表。质量输入 (genomeInfo) 从 Phase 8 的
    // QC 表构建, 无外部数据库依赖。产出:
    //   representatives = tuple(meta, mag_id, mag_fasta)  (形状同输入, Phase 10 消费)
    //   membership      = mag_id \t sample \t representative_mag_id (无表头)
    //   catalog/clusters = 代表 MAG FASTA 目录 / dRep 原始数据表 (Phase 11/13 消费)
    //
    // --skip_dereplication 与 --skip_mag_qc 交互守卫在 DREP 子工作流内;
    // skip 时代表集 = 全部 qualified MAG, 成员表为恒等映射 (等价全量分类)。
    // ---------------------------------------------------------------------
    DREP(ch_qualified_mags, ch_qc_table)
    ch_representatives = DREP.out.representatives
    ch_membership      = DREP.out.membership

    // ---------------------------------------------------------------------
    // Phase 10: MAG 分类 (GTDB-Tk classify_wf → mag_taxonomy.tsv)
    //
    // 执行顺序为**先去冗余后分类**: 本阶段只分类 Phase 9 的代表 MAG, 冗余
    // 成员 MAG 的分类由 dRep 成员表回填 (见 TAXONOMY 子工作流)。
    //
    // --skip_taxonomy 与 --gtdbtk_db 守卫在 TAXONOMY 子工作流内。
    // ---------------------------------------------------------------------
    TAXONOMY(ch_representatives, ch_membership)
    ch_taxonomy_table = TAXONOMY.out.taxonomy_table
    ch_multiqc_files  = ch_multiqc_files.mix(TAXONOMY.out.multiqc_files)  // Phase 15: mag_taxonomy.tsv

    // ---------------------------------------------------------------------
    // Phase 11: 基因预测 (Prodigal → genes.fna / proteins.faa / GFF)
    //
    // 逐代表 MAG 各调用一次 (-p meta 宏基因组模式), 无外部数据库依赖。
    // 上游 representatives 通道本身已是 tuple(meta, mag_id, mag_fasta),
    // 直接映射, 无需集合级聚合。proteins 通道是 Phase 12 三支功能注释
    // (DIAMOND / eggNOG-mapper / RGI) 的唯一输入。
    // --skip_gene_prediction 守卫在 GENE_PREDICTION 子工作流内。
    // ---------------------------------------------------------------------
    GENE_PREDICTION(ch_representatives)
    ch_proteins = GENE_PREDICTION.out.proteins
    ch_genes    = GENE_PREDICTION.out.genes
    ch_gff      = GENE_PREDICTION.out.gff

    // ---------------------------------------------------------------------
    // Phase 12: 功能注释 (DIAMOND blastp ∥ eggNOG-mapper ∥ RGI)
    //
    // 三支彼此独立的并行分支, 唯一输入是 Phase 11 的 proteins 通道
    // (tuple(meta, mag_id, proteins.faa), 每个代表 MAG 一条), 逐 MAG 直接
    // 映射, 无集合级聚合。每支产出原始工具输出 + 一张 (meta_id, mag_id,
    // gene) 键控的解析汇总表 —— Phase 14 总表 join 的键。
    // skip/数据库守卫在 ANNOTATION 子工作流内 (Phase 4/8/10/11 模式)。
    // ---------------------------------------------------------------------
    ANNOTATION(ch_proteins)
    ch_diamond_table = ANNOTATION.out.diamond_table
    ch_eggnog_table  = ANNOTATION.out.eggnog_table
    ch_rgi_table     = ANNOTATION.out.rgi_table

    // ---------------------------------------------------------------------
    // Phase 13: MAG 丰度 (CoverM genome → 跨样本丰度矩阵)
    //
    // 集合级批处理: 全部代表 MAG FASTA (Phase 9) + 全部样本的 Phase 6 排序
    // BAM 一次调用产出矩阵 (行 = mag_id, 列 = 样本)。BAM 直接复用 Phase 6
    // 产物 —— MAG FASTA 的 contig header 原样保留, 与 BAM 的 contig 名一致,
    // 无需重新比对。无外部数据库依赖 (CoverM 自带比对器)。
    //
    // --skip_mapping / --skip_assembly 时没有 BAM, --skip_binning 时没有代表
    // MAG, 均在此明确告警并跳过 (与 Phase 5/6/7 的 skip 交互守卫同模式);
    // 其余上游组合导致通道为空时, 子工作流内聚合不发射, 正常收尾。
    // ---------------------------------------------------------------------
    ch_abundance = Channel.empty()

    if (!params.skip_abundance) {
        if (params.skip_mapping) {
            log.warn "MAG 丰度已跳过: --skip_mapping 时无排序 BAM 可用, CoverM 无从运行 (不产出 13_abundance)。"
        }
        else if (params.skip_assembly) {
            log.warn "MAG 丰度已跳过: --skip_assembly 时无 contigs, 也没有可复用的比对 BAM (不产出 13_abundance)。"
        }
        else if (params.skip_binning) {
            log.warn "MAG 丰度已跳过: --skip_binning 时无 MAG 可计算丰度 (不产出 13_abundance)。"
        }
        else {
            ABUNDANCE(MAPPING.out.bam, ch_representatives)
            ch_abundance     = ABUNDANCE.out.abundance
            ch_multiqc_files = ch_multiqc_files.mix(ABUNDANCE.out.multiqc_files)  // Phase 15: mag_abundance.tsv
        }
    }
    else {
        log.warn "MAG 丰度已跳过 (--skip_abundance), 不产出 13_abundance 丰度矩阵。"
    }

    // ---------------------------------------------------------------------
    // Phase 14: 结果整合 (各 Phase 汇总表 → 14_integrated 核心表)
    //
    // 只消费上游通道, 不重做任何计算: QC/分类/成员表/注释/丰度/代表 MAG/
    // 组装汇总 join 成两张核心表 (mag_metadata.tsv / mag_functional_annotation.tsv)
    // + 引用拷贝 (成员表/组装汇总)。每个输入表都可能因对应 --skip_* 为空
    // 通道 —— INTEGRATION 子工作流内以 0 字节哨兵兜底 (ifEmpty 只接受具体
    // 值, 详见文档), 缺失表的对应列留空, 行集合由存在表的并集
    // 决定。pathway_db 为可选外部依赖 (KO→pathway 映射), 未提供时 Pathway
    // 列留空 (不虚构)。上游组合导致代表 MAG 通道为空时 (如 --skip_binning),
    // 子工作流内聚合不发射, 本阶段不调度任务 (正常收尾)。
    // ---------------------------------------------------------------------
    if (!params.skip_integration) {
        INTEGRATION(
            ch_representatives, ch_membership,
            ch_qc_table, ch_taxonomy_table, ch_abundance,
            ch_diamond_table, ch_eggnog_table, ch_rgi_table,
            ch_assembly_summary
        )

        // 两张核心表供 Phase 15 MultiQC 消费 (混入 ch_multiqc_files)
        ch_multiqc_files = ch_multiqc_files.mix(INTEGRATION.out.multiqc_files)
    }
    else {
        log.warn "结果整合已跳过 (--skip_integration), 不产出 14_integrated 汇总表。"
    }

    // ---------------------------------------------------------------------
    // Phase 21: 结果可视化 (画图消费层 —— 只读各 Phase 已产出 TSV, 不重做计算)
    //
    // 每个 plot process 以其输入通道是否为空独立调度: 对应 Phase 被 skip 时
    // 该图不调度、figures/ 不建。database-dependent 图 (除 MAG 丰度热图外)
    // 真实数据本地无库, 代码由 stub/合成表验证"能画", 真实图待库 。
    // ---------------------------------------------------------------------
    PLOTTING(
        ch_bin_summary, ch_qc_table, ch_taxonomy_table, ch_membership,
        ch_abundance, ch_bracken_merged, ch_beta_diversity,
        ch_pathabundance_merged
    )

    // ---------------------------------------------------------------------
    // 汇总 QC 报告 (Phase 15)
    //
    // ch_multiqc_files 汇集各阶段 QC 产物 (Phase 3-8 / 10 / 13 / 14)。
    // Phase 9 (dRep) / 11 (Prodigal) / 12 (原始注释表) 无 MultiQC 可解析
    // 产物, 注释三表经 Phase 14 聚合的 mag_functional_annotation.tsv 进入
    // 报告。自定义 TSV 解析配置经 ch_multiqc_config 传入 MULTIQC。
    // ---------------------------------------------------------------------
    ch_multiqc_config = file("${projectDir}/modules/local/qc/multiqc_config.yaml", checkIfExists: true)

    if (!params.skip_multiqc) {
        MULTIQC(ch_multiqc_files.collect(), ch_multiqc_config)
    }
}
