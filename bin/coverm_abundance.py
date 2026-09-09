#!/usr/bin/env python3
"""
coverm_abundance.py — Phase 13 COVERM 封装: BAM 过滤 + 集合级调用 + 矩阵规范化

职责 (三合一, 全部数据留在进程任务目录内):
  1. BAM 过滤: 按 contig 名交集筛掉与任何 MAG 都无重叠的 BAM。CoverM 对
     「与任何 genome 都无 contig 交集」的 BAM 直接报错 ("There are no found
     reference sequences that are a part of a genome"), 而不是产出全零列 ——
     0-bin 样本 (该样本的 BAM 里没有任何 contig 属于某个 MAG) 必然触发。
     无交集的 BAM 不进 --bam-files, 其样本列由规范化步骤补 0。
  2. 调用 CoverM genome: 集合级**单次调用** (全部有交集 BAM + 全部代表 MAG),
     直接产出原始矩阵。
  3. 矩阵规范化 (Phase 14 join 的前提):
       列 = bam manifest 中的列名 (单组装器 = 样本 ID; --assembler both 时
            同一样本有两套 BAM, 子工作流已消歧为 <sample>.<assembler>),
            被剔除 BAM 的样本列补 0, 保证列集合 = 全部样本
       行 = mag_id (MAG FASTA 去 .fa 后缀), 按名排序保证 -resume 哈希稳定;
            剔除 unmapped 伪行 (未落在任何代表 MAG 上的 reads 占比)
       数值 = 表头带 "(%)" 的列 (相对丰度类方法) 由 0-100 百分比换算为
            0-1 小数, 其余方法 (mean / rpkm / tpm 等) 原样保留

─── BAM 表头读取 (不引入 samtools 依赖) ───────────────────────────────────
    BAM = BGZF (gzip 兼容), 流结构: 魔数 "BAM\\x01"(4B) + l_text(4B) + 文本头
    (l_text 字节)。用 gzip + struct 按字节偏移精确截取 —— 不能按行正则解析
    (l_text 的 4 个字节可能恰含换行)。@SQ 的 SN 即参考序列名, 与 CoverM 的
    匹配口径一致 (CoverM 默认只比较 contig 名空白前的部分, MAG FASTA 侧同理
    取 header 第一段)。

─── 参数 ──────────────────────────────────────────────────────────────────
    --bam-manifest  <BAM basename> \\t <列名> 清单 (collectFile 物化, 排序)
    --genomes       全部代表 MAG FASTA (<mag_id>.fa)
    --method        CoverM --methods 值
    --threads       CoverM -t
    --coverm-args   追加给 coverm genome 的额外参数 (shlex 拆分)
    --output / --raw-output  规范矩阵 / CoverM 原始矩阵

仅依赖标准库 (gzip/struct/subprocess/shlex/csv/argparse), 便于在任何环境
(conda / 容器) 与单元测试中复用。
"""

import argparse
import csv
import gzip
import shlex
import struct
import subprocess
import sys
from pathlib import Path

BAM_MAGIC = b"BAM\x01"


def bam_reference_names(bam_path: Path) -> set:
    """BAM 文本头中 @SQ 的 SN 集合 (参考序列名)。"""
    with gzip.open(bam_path, "rb") as fh:
        if fh.read(4) != BAM_MAGIC:
            raise ValueError(f"{bam_path}: 不是有效的 BAM 文件 (魔数不匹配)")
        (l_text,) = struct.unpack("<I", fh.read(4))
        header_text = fh.read(l_text).decode("ascii", errors="replace")

    names = set()
    for line in header_text.splitlines():
        if line.startswith("@SQ"):
            for field in line.split("\t"):
                if field.startswith("SN:"):
                    names.add(field[3:])
                    break
    return names


def mag_contig_names(fasta_files: list) -> set:
    """MAG FASTA header 第一段集合 (与 CoverM 匹配口径一致)。"""
    names = set()
    for path in fasta_files:
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    names.add(line[1:].split()[0])
    return names


def normalize(raw_path: Path, out_path: Path, bam2col: dict, col_order: list):
    """CoverM 原始矩阵 → 键控规范化矩阵 (行 = mag_id, 列 = 全部样本)。"""
    with open(raw_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # 列名: "<BAM名去.bam> <方法显示名> (%)" → 映射为样本列名
        col_map, pct = {}, []
        for i, cell in enumerate(header[1:], start=1):
            base = cell.split(" ", 1)[0]
            if base.endswith(".bam"):
                base = base[:-4]
            pct.append(cell.endswith("(%)"))
            col_map[i] = bam2col.get(base, base)

        rows = []
        for line in fh:
            cells = line.rstrip("\n").split("\t")
            name = cells[0]
            if name.endswith(".fa"):
                name = name[:-3]
            if name == "unmapped":
                continue
            values = {
                col_map[i]: (str(float(v) / 100) if pct[i - 1] else v)
                for i, v in enumerate(cells[1:], start=1)
            }
            rows.append((name, values))

    # 按 mag_id 排序: 输出内容确定, -resume 哈希稳定
    rows.sort(key=lambda r: r[0])

    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["MAG_ID", *col_order])
        for name, values in rows:
            writer.writerow([name, *[values.get(c, "0") for c in col_order]])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Phase 13 CoverM 封装")
    ap.add_argument("--bam-manifest", required=True, type=Path)
    ap.add_argument("--genomes", nargs="+", required=True, type=Path)
    ap.add_argument("--method", required=True)
    ap.add_argument("--threads", type=int, required=True)
    ap.add_argument("--coverm-args", default="")
    ap.add_argument("--output", default=Path("mag_abundance.tsv"), type=Path)
    ap.add_argument("--raw-output", default=Path("mag_abundance_raw.tsv"), type=Path)
    args = ap.parse_args(argv)

    # 1. BAM 清单 → (basename, 列名), 保持清单顺序即列序
    # 键统一为「basename 去 .bam」: CoverM 表头列名即此口径 (实测 0.8.0)
    manifest = []
    with open(args.bam_manifest) as fh:
        for line in fh:
            bam_name, col = line.rstrip("\n").split("\t")
            manifest.append((bam_name, col))
    bam2col = {
        (b[:-4] if b.endswith(".bam") else b): col
        for b, col in manifest
    }
    col_order = [col for _, col in manifest]

    # 2. 按 contig 名交集过滤 BAM
    mag_names = mag_contig_names(args.genomes)
    included = [(b, c) for b, c in manifest if bam_reference_names(Path(b)) & mag_names]
    if not included:
        print("ERROR: 没有任何 BAM 与代表 MAG 存在 contig 交集, CoverM 无法运行。",
              file=sys.stderr)
        return 1

    # 3. 集合级单次调用 CoverM genome
    cmd = [
        "coverm", "genome",
        "--bam-files", *[b for b, _ in included],
        "--genome-fasta-files", *[str(g) for g in args.genomes],
        "--methods", args.method,
        "--output-file", str(args.raw_output),
        "--threads", str(args.threads),
    ]
    if args.coverm_args.strip():
        cmd.extend(shlex.split(args.coverm_args))
    subprocess.run(cmd, check=True)

    # 4. 键控规范化
    normalize(args.raw_output, args.output, bam2col, col_order)
    return 0


if __name__ == "__main__":
    sys.exit(main())
