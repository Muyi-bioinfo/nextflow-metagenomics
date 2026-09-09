#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""integrate_metadata.py — Phase 14 整合: 各 Phase 汇总表 → mag_metadata.tsv

输入:
  --manifest    三列 TSV (meta_id, mag_id, fasta_path), 无表头, 由 INTEGRATION
                子工作流 collectFile 产出 (sort: true, -resume 哈希稳定) ——
                代表 MAG 的权威清单, 输出行集合与之严格一致 (即 dRep 代表集)。
                fasta_path 为上游 work 目录绝对路径, resolve() 双路径打开:
                真实运行的暂存 basename 优先, -stub-run 不暂存时回退 manifest
                绝对路径 (Phase 12 模式)。
  --membership  三列 TSV (mag_id, sample, representative_mag_id), 无表头,
                Phase 9 成员表 (dRep 聚类或恒等映射)。Sample 列取自此表
                (与 Phase 10 分类回填同源); 清单 MAG 必须是成员表中的代表
                行 (rep == mag_id), 否则报错。
  --qc          可选: Phase 8 mag_qc.tsv (表头 meta_id/mag_id/Name/
                Completeness/Contamination)。缺失或 0 字节 = CheckM2 未运行,
                Completeness/Contamination 列留空。QC 表含全部输入 MAG (含
                被阈值过滤的), 只取清单中的代表 MAG 行, 多余行属预期不报错。
  --taxonomy    可选: Phase 10 mag_taxonomy.tsv (表头 sample/mag_id/
                rep_mag_id/domain/phylum/class/order/family/genus/species)。
                取代表 MAG 自身行 (rep_mag_id 为空); GTDB_taxonomy 由各 rank
                按 d__;p__;c__;o__;f__;g__;s__ 前缀重建为完整分类串 (缺级
                rank 省略, 全缺为空串), Species 为 s__ 一级 (不含前缀)。
  --abundance   可选: Phase 13 mag_abundance.tsv (表头 MAG_ID + 样本 ID 列,
                数值 0-1 字符串)。每样本一列宽表化并入输出 (列名 = 样本 ID,
                顺序同矩阵列序), 代表 MAG 无行或单元格缺失时补 0。

Genome_size / GC 由各代表 MAG FASTA 现场计算 (Python 标准库, 不引入新工具):
Genome_size = 序列总长 (bp, 整数), GC = (G+C)/总长百分比 (保留两位小数)。

输出 mag_metadata.tsv (列序固定):
  MAG_ID  Sample  Completeness  Contamination  Genome_size  GC
  GTDB_taxonomy  Species  [<样本 ID> ...]   (丰度表存在时追加样本列)

行排序 = manifest 顺序 (collectFile sort:true 已按内容排序, 输出确定)。
"""
import argparse
import csv
import os
import sys

BASE_COLUMNS = ["MAG_ID", "Sample", "Completeness", "Contamination",
                "Genome_size", "GC", "GTDB_taxonomy", "Species"]

# 分类 rank: (GTDB-Tk 前缀, 表列名) —— 与 taxonomy_summary.py 一致
RANKS = (
    ("d__", "domain"),
    ("p__", "phylum"),
    ("c__", "class"),
    ("o__", "order"),
    ("f__", "family"),
    ("g__", "genus"),
    ("s__", "species"),
)
RANK_COLUMNS = [name for _, name in RANKS]


def resolve(path):
    """解析 manifest 中的文件路径: 任务目录暂存名 (basename) 优先, 回退
    绝对路径 (Phase 12 模式: 真实运行暂存 basename 命中; -stub-run 未被
    stub 引用的输入不暂存, 回退 manifest 绝对路径, 上游 work 文件仍在)。
    两者都找不到才报错。"""
    base = os.path.basename(path)
    if os.path.exists(base):
        return base
    if os.path.exists(path):
        return path
    raise SystemExit(f"ERROR: cannot resolve {path!r}: neither staged basename "
                     f"nor manifest absolute path exists")


def load_manifest(path):
    """读代表 MAG 清单 (meta_id \\t mag_id \\t fasta_path, 无表头)。"""
    rows = []
    seen = set()
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            fields = line.rstrip("\n").split("\t")
            if not fields or not fields[0]:
                continue
            if len(fields) < 3:
                raise SystemExit(f"manifest {path} line {lineno} malformed: {line!r}")
            meta_id, mag_id, fasta_path = fields[0], fields[1], fields[2]
            if mag_id in seen:
                raise SystemExit(f"duplicate mag_id in manifest: {mag_id}")
            seen.add(mag_id)
            rows.append((meta_id, mag_id, fasta_path))
    return rows


def load_membership(path):
    """读成员表 (mag_id \\t sample \\t representative_mag_id, 无表头)。"""
    membership = {}
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            fields = line.rstrip("\n").split("\t")
            if not fields or not fields[0]:
                continue
            if len(fields) < 3:
                raise SystemExit(f"membership {path} line {lineno} malformed: {line!r}")
            mag_id, sample, rep = fields[0], fields[1], fields[2]
            if mag_id in membership:
                raise SystemExit(f"duplicate mag_id in membership: {mag_id}")
            membership[mag_id] = (sample, rep)
    return membership


def load_qc(path):
    """读 QC 表 (可选)。缺失/0 字节 → {}; 按 mag_id 键控, 多余 MAG (被
    阈值过滤的) 保留在返回结果中, 由调用方按清单取用。"""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        print(f"note: qc table {path or '(none)'} absent or empty — QC columns blank",
              file=sys.stderr)
        return {}
    qc = {}
    with open(path, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for required in ("mag_id", "Completeness", "Contamination"):
            if required not in reader.fieldnames:
                raise SystemExit(f"qc table {path} lacks required column {required!r}")
        for row in reader:
            mag_id = row["mag_id"].strip()
            if mag_id in qc:
                raise SystemExit(f"duplicate mag_id in qc table: {mag_id}")
            qc[mag_id] = (row.get("Completeness", "").strip(),
                          row.get("Contamination", "").strip(),
                          row.get("meta_id", "").strip())
    return qc


def load_taxonomy(path):
    """读分类表 (可选)。缺失/0 字节 → None (表不存在, 列留空); 表存在时
    取代表 MAG 自身行 (rep_mag_id 为空), 键 = mag_id, 值 = (GTDB_taxonomy
    串, species) —— 两者必须区分: 表存在而代表 MAG 无行属接线断裂 (Phase
    10 契约: 分类失败也保留空分类行), 表不存在则整列留空。"""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        print(f"note: taxonomy table {path or '(none)'} absent or empty — "
              f"taxonomy columns blank", file=sys.stderr)
        return None
    taxonomy = {}
    with open(path, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for required in ["mag_id", "rep_mag_id"] + RANK_COLUMNS:
            if required not in reader.fieldnames:
                raise SystemExit(f"taxonomy table {path} lacks required column {required!r}")
        for row in reader:
            mag_id = row["mag_id"].strip()
            rep = row["rep_mag_id"].strip()
            if rep:
                continue  # 冗余成员回填行 —— 元数据表行集合 = 代表 MAG
            tokens = [f"{prefix}{row[name].strip()}"
                      for prefix, name in RANKS
                      if row[name].strip()]
            gtdb_string = ";".join(tokens)
            if mag_id in taxonomy:
                raise SystemExit(f"duplicate representative row in taxonomy table: {mag_id}")
            taxonomy[mag_id] = (gtdb_string, row["species"].strip())
    return taxonomy


def load_abundance(path):
    """读丰度矩阵 (可选)。缺失/0 字节 → (None, None); 否则返回
    (样本列名列表, {mag_id: {样本: 值}}), 缺失值由调用方补 0。"""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        print(f"note: abundance matrix {path or '(none)'} absent or empty — "
              f"abundance columns omitted", file=sys.stderr)
        return None, None
    abundance = {}
    with open(path, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or reader.fieldnames[0] != "MAG_ID":
            raise SystemExit(f"abundance matrix {path} must have MAG_ID as first column")
        samples = list(reader.fieldnames[1:])
        for row in reader:
            mag_id = row["MAG_ID"].strip()
            if mag_id in abundance:
                raise SystemExit(f"duplicate MAG_ID in abundance matrix: {mag_id}")
            # 空单元格与缺失列一律补 0 (矩阵形状 = 行 mag_id × 列 全部样本)
            abundance[mag_id] = {s: row.get(s, "").strip() or "0" for s in samples}
    return samples, abundance


def fasta_stats(path):
    """计算 FASTA 的总长 (bp) 与 GC 百分比 (两位小数字符串)。"""
    size, gc = 0, 0
    with open(resolve(path), encoding="ascii", errors="strict") as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq = line.strip()
            size += len(seq)
            gc += seq.count("G") + seq.count("g") + seq.count("C") + seq.count("c")
    if size == 0:
        raise SystemExit(f"ERROR: FASTA {path} has no sequence — cannot compute "
                         f"Genome_size/GC")
    return str(size), f"{gc * 100.0 / size:.2f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--membership", required=True)
    ap.add_argument("--qc", default="")
    ap.add_argument("--taxonomy", default="")
    ap.add_argument("--abundance", default="")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    manifest = load_manifest(args.manifest)
    membership = load_membership(args.membership)

    # 清单 MAG 必须是成员表中的代表行 —— 行集合 = dRep 代表集的键一致性检查
    for meta_id, mag_id, _fasta in manifest:
        entry = membership.get(mag_id)
        if entry is None:
            raise SystemExit(f"representative {mag_id!r} missing from membership — "
                             f"wiring broken")
        sample, rep = entry
        if rep != mag_id:
            raise SystemExit(f"{mag_id!r} is listed as a non-representative "
                             f"(rep={rep!r}) in membership — wiring broken")

    qc = load_qc(args.qc)
    taxonomy = load_taxonomy(args.taxonomy)
    samples, abundance = load_abundance(args.abundance)

    columns = list(BASE_COLUMNS) + (samples or [])

    rows = []
    for meta_id, mag_id, fasta_path in manifest:
        sample, _rep = membership[mag_id]
        completeness, contamination, qc_meta_id = qc.get(mag_id, ("", "", ""))
        if qc_meta_id and qc_meta_id != meta_id:
            raise SystemExit(f"meta_id mismatch for {mag_id!r}: manifest "
                             f"{meta_id!r} vs qc table {qc_meta_id!r} — wiring broken")
        if taxonomy is None:
            gtdb_string, species = "", ""
        else:
            # 分类表存在时每个代表 MAG 都必须有自身行 (Phase 10 契约: 分类
            # 失败也保留空分类行), 缺失即接线断裂
            if mag_id not in taxonomy:
                raise SystemExit(f"representative {mag_id!r} missing from taxonomy "
                                 f"table — wiring broken")
            gtdb_string, species = taxonomy[mag_id]
        genome_size, gc = fasta_stats(fasta_path)
        row = [mag_id, sample, completeness, contamination, genome_size, gc,
               gtdb_string, species]
        if samples:
            values = abundance.get(mag_id, {})
            row.extend(values.get(s, "0") for s in samples)
        rows.append(row)

    with open(args.output, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


if __name__ == "__main__":
    main()
