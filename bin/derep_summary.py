#!/usr/bin/env python3
"""derep_summary.py — dRep 聚类输出 → 成员表 (Phase 9)

输入:
  --manifest      genome_manifest.tsv:  mag_id \\t sample (全部 qualified MAG,
                  由 DREP 子工作流 collectFile 生成, 无表头)
  --cdb           dRep data_tables/Cdb.csv (genome,secondary_cluster)
  --dereplicated  dRep dereplicated_genomes/ 目录 (每簇一个代表 FASTA)

输出:
  mag_membership.tsv   mag_id \\t sample \\t representative_mag_id (无表头)

mag_id 与 dRep genome 名的对应: MAG FASTA 以 <mag_id>.fa 命名 (Phase 7 稳定
MAG ID 保证), dRep 输出表的 genome 列是去扩展名的 basename, 即 mag_id 本身。

代表判定: dereplicated_genomes/ 中的 FASTA 即 dRep 选出的簇代表; Cdb.csv 的
secondary_cluster 给出每个 MAG 所属的簇, 簇内恰有一个代表。任何映射不一致
(重复 / 缺失 / 一簇多代表) 都直接报错, 不静默回退 —— 成员表是 Phase 10
分类回填的依据, 错一行即错一片。
"""
import argparse
import csv
import sys
from pathlib import Path

FASTA_EXTS = (".fa", ".fasta", ".fna")


def strip_fasta_ext(name: str) -> str:
    for ext in FASTA_EXTS:
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--cdb", required=True, type=Path)
    ap.add_argument("--dereplicated", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    # manifest: mag_id → sample (保持插入序, 重复 mag_id 直接报错)
    manifest = {}
    with args.manifest.open() as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                sys.exit(f"ERROR: manifest line {lineno} has < 2 columns: {line!r}")
            mag_id, sample = cols[0].strip(), cols[1].strip()
            if not mag_id or not sample:
                sys.exit(f"ERROR: manifest line {lineno}: empty mag_id or sample: {line!r}")
            if mag_id in manifest:
                sys.exit(f"ERROR: duplicate mag_id in manifest: {mag_id}")
            manifest[mag_id] = sample

    # Cdb.csv: genome → secondary_cluster (按位置取列, 不依赖列名)
    clusters = {}
    with args.cdb.open() as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if header is None:
            sys.exit("ERROR: Cdb.csv is empty (no header)")
        for lineno, row in enumerate(reader, 2):
            if not row or all(c.strip() == "" for c in row):
                continue
            if len(row) < 2:
                sys.exit(f"ERROR: Cdb.csv line {lineno} has < 2 columns")
            genome = row[0].strip()
            if genome in clusters:
                sys.exit(f"ERROR: duplicate genome in Cdb.csv: {genome}")
            clusters[genome] = row[1].strip()

    # 代表: dereplicated_genomes/ 中的 FASTA 文件
    if not args.dereplicated.is_dir():
        sys.exit(f"ERROR: dereplicated_genomes dir not found: {args.dereplicated}")
    winners = set()
    for p in sorted(args.dereplicated.iterdir()):
        if p.is_file() and p.name.endswith(FASTA_EXTS):
            winners.add(strip_fasta_ext(p.name))
    if not winners:
        sys.exit("ERROR: no representative FASTA in dereplicated_genomes/")

    # 每个簇的代表必须唯一; 代表不在 Cdb 时视为单簇 (防御性容错,
    # dRep 正常输出不会触发 —— 正常输出中代表一定在 Cdb 里)
    cluster_rep = {}
    for w in winners:
        if w not in manifest:
            sys.exit(f"ERROR: representative {w} not in manifest (genome name mismatch)")
        c = clusters.get(w)
        if c is not None:
            if cluster_rep.get(c, w) != w:
                sys.exit(f"ERROR: secondary_cluster {c} has multiple representatives: "
                         f"{cluster_rep[c]} and {w}")
            cluster_rep[c] = w

    rows = []
    for mag_id, sample in manifest.items():
        if mag_id in winners:
            rep = mag_id
        else:
            c = clusters.get(mag_id)
            if c is None:
                sys.exit(f"ERROR: genome {mag_id} missing from Cdb.csv — cannot resolve representative")
            rep = cluster_rep.get(c)
            if rep is None:
                sys.exit(f"ERROR: secondary_cluster {c} (member {mag_id}) has no "
                         f"representative in dereplicated_genomes/")
        rows.append((mag_id, sample, rep))

    rows.sort()
    with args.output.open("w") as fh:
        for mag_id, sample, rep in rows:
            fh.write(f"{mag_id}\t{sample}\t{rep}\n")


if __name__ == "__main__":
    main()
