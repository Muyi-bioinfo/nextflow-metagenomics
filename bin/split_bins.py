#!/usr/bin/env python3
"""
split_bins.py — 把分箱器输出的 bin 目录拆分成独立命名的 MAG FASTA

用途 (Phase 7):
    MetaBAT2 以「输出前缀」的方式产出 <prefix>.1.fa / <prefix>.2.fa / ...
    这些文件名只在单个组装单元内唯一 —— 一旦多个样本、多个组装器、多个分箱器
    的结果汇到一起 (Phase 8 CheckM2 / Phase 10 dRep), 就会互相撞名。
    本脚本把每个 bin 复制成带全局唯一 MAG ID 的文件, 并输出一张 per-MAG 统计表。

─── MAG ID 的构成与稳定性 ─────────────────────────────────────────────────
    mag_id = <组装单元 id>.<assembler>.<binner>.<bin 序号, 三位补零>
    例:      S01.megahit.metabat2.001

    四个字段各自解决一类撞名:
      组装单元 id   不同样本 (或 coassembly 组) 的 bin
      assembler     --assembler both 时同一样本有两套独立 contigs
      binner        V2 接入 MaxBin2 / CONCOCT / DAS Tool 后同一套 contigs 有多份分箱
      bin 序号      同一次分箱内的不同 bin

    序号**沿用分箱器自己给出的编号**, 不重新连续编号。原因:
      1. 可追溯 —— MAG ID 能直接对回 07_binning/metabat2/ 下的原始文件;
      2. 稳定 —— 重新编号会让「某个 bin 因阈值变化消失」连带改变它之后所有
         MAG 的 ID; 沿用原编号则只影响消失的那一个。

    这套 ID 的稳定性还依赖分箱器本身可重复。MetaBAT2 需固定 --seed
    (见 modules/local/binning/metabat2.nf 与 params.metabat2_seed), 否则
    同样的输入可能产出不同的 bin 划分, MAG ID 也就跟着变。

─── 不改动序列, 也不改 contig header ───────────────────────────────────────
    MAG FASTA 里的 contig header 原样保留。改写 header 会切断 MAG 与 05_assembly/
    contigs、06_mapping/ 深度矩阵之间的对应关系 —— Phase 13 (CoverM MAG 丰度)
    正是靠 contig 名把 MAG 映射回比对结果的。

─── 哪些文件算 bin ────────────────────────────────────────────────────────
    只接受文件名倒数第二段是**整数**的文件 (<prefix>.<N>.fa)。
    这样 MetaBAT2 --unbinned 产出的 <prefix>.unbinned.fa, 以及
    <prefix>.tooShort.fa / <prefix>.lowDepth.fa 之类的非 bin 产物会被自动排除,
    不需要维护一张后缀黑名单。

─── 0 个 bin 是合法结果 ───────────────────────────────────────────────────
    分箱不出任何 bin 时本脚本正常退出 (退出码 0), 只写出一张仅含表头的统计表,
    并在 stderr 说明。这不是失败 —— 低复杂度或低深度数据下分箱器形成不了满足
    最小 bin 尺寸的簇是真实结果。伪造一个空 MAG 才是错的。
"""

import argparse
import shutil
import sys
from pathlib import Path

# 统计表列顺序 (Phase 14 汇总与 BIN_SUMMARY 拼接均依赖此顺序)
COLUMNS = [
    "mag_id",
    "assembly_unit",
    "assembler",
    "assembly_mode",
    "binner",
    "source_bin",
    "n_contigs",
    "total_bp",
    "largest_contig_bp",
    "mean_contig_bp",
    "gc_percent",
    "n_bases",
]

# 认作 FASTA 的扩展名 (分箱器实际只产 .fa, 另两个为容错)
FASTA_SUFFIXES = {".fa", ".fna", ".fasta"}


def parse_args():
    p = argparse.ArgumentParser(description="把分箱器的 bin 目录拆分为带稳定 MAG ID 的独立 FASTA")
    p.add_argument("--bins-dir", required=True, help="分箱器输出目录 (内含 <prefix>.<N>.fa)")
    p.add_argument("--unit-id", required=True, help="组装单元 id (meta.id)")
    p.add_argument("--assembler", required=True, help="产出 contigs 的组装器 (megahit | metaspades)")
    p.add_argument("--assembly-mode", required=True, help="single | coassembly")
    p.add_argument("--binner", required=True, help="分箱器名 (metabat2 | ...)")
    p.add_argument("--outdir", required=True, help="MAG FASTA 输出目录")
    p.add_argument("--summary", required=True, help="per-MAG 统计表输出路径 (TSV)")
    return p.parse_args()


def bin_index(path: Path):
    """
    从 <prefix>.<N>.fa 提取 bin 序号 N。

    返回 int, 或 None 表示这不是一个 bin 文件 (倒数第二段不是整数)。
    """
    if path.suffix.lower() not in FASTA_SUFFIXES:
        return None

    parts = path.name[: -len(path.suffix)].split(".")
    if len(parts) < 2:
        return None

    try:
        return int(parts[-1])
    except ValueError:
        return None


def read_fasta_lengths(path: Path):
    """
    单遍扫描一个 FASTA, 返回 (每条 contig 的长度列表, GC 计数, N 计数)。

    GC 只在 ACGT 之中统计: 把 N 计入分母会让高 N 的 MAG 的 GC 被系统性低估。
    """
    lengths = []
    gc = 0
    n_bases = 0
    acgt = 0
    current = 0
    seen_header = False

    with path.open("r") as fh:
        for line in fh:
            if line.startswith(">"):
                seen_header = True
                if current:
                    lengths.append(current)
                current = 0
                continue
            seq = line.strip()
            current += len(seq)
            for base in seq.upper():
                if base in ("G", "C"):
                    gc += 1
                    acgt += 1
                elif base in ("A", "T"):
                    acgt += 1
                elif base == "N":
                    n_bases += 1

    if current:
        lengths.append(current)

    if not seen_header:
        raise ValueError(f"{path.name} 不是 FASTA (没有任何以 '>' 开头的行)")

    return lengths, gc, acgt, n_bases


def main():
    args = parse_args()

    bins_dir = Path(args.bins_dir)
    if not bins_dir.is_dir():
        sys.exit(f"ERROR: bin 目录不存在: {bins_dir}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 按分箱器给出的序号排序, 保证输出顺序与 MAG ID 顺序一致且可重复
    candidates = []
    for path in sorted(bins_dir.iterdir()):
        if not path.is_file():
            continue
        idx = bin_index(path)
        if idx is None:
            print(f"跳过非 bin 文件: {path.name}", file=sys.stderr)
            continue
        candidates.append((idx, path))
    candidates.sort(key=lambda x: x[0])

    rows = []
    for idx, path in candidates:
        mag_id = f"{args.unit_id}.{args.assembler}.{args.binner}.{idx:03d}"

        lengths, gc, acgt, n_bases = read_fasta_lengths(path)
        if not lengths:
            # bin 文件存在但没有序列: 分箱器不该产出这种文件, 不静默放过
            sys.exit(f"ERROR: {path.name} 内没有任何序列, 无法作为 MAG 输出。")

        total = sum(lengths)
        # 序列内容原样复制 —— header 保持不变 (见文件头说明)
        shutil.copyfile(path, outdir / f"{mag_id}.fa")

        rows.append(
            {
                "mag_id": mag_id,
                "assembly_unit": args.unit_id,
                "assembler": args.assembler,
                "assembly_mode": args.assembly_mode,
                "binner": args.binner,
                "source_bin": path.name,
                "n_contigs": len(lengths),
                "total_bp": total,
                "largest_contig_bp": max(lengths),
                "mean_contig_bp": round(total / len(lengths), 1),
                "gc_percent": round(100.0 * gc / acgt, 2) if acgt else "NA",
                "n_bases": n_bases,
            }
        )

    with Path(args.summary).open("w") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            out.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")

    if not rows:
        # 合法结果, 不是错误 —— 见文件头「0 个 bin 是合法结果」
        print(
            f"未在 {bins_dir} 找到任何 bin, {args.summary} 只含表头。"
            f" 若这不符合预期, 检查 {args.binner} 的最小 bin 尺寸与最小 contig 长度阈值。",
            file=sys.stderr,
        )
    else:
        print(f"拆分出 {len(rows)} 个 MAG -> {outdir}/", file=sys.stderr)
        for row in rows:
            print(
                f"  {row['mag_id']}\t{row['n_contigs']} contigs\t{row['total_bp']} bp"
                f"\tGC {row['gc_percent']}%",
                file=sys.stderr,
            )


if __name__ == "__main__":
    main()
