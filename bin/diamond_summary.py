#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Parse per-MAG DIAMOND blastp outputs into one keyed best-hit table.

输入:
  --manifest  四列 TSV (meta_id, mag_id, diamond_tsv, proteins_faa), 无表头,
              由 ANNOTATION 子工作流 collectFile 产出 —— 本批次全部注释 MAG
              的权威清单。文件路径经 resolve() 解析: 真实运行的暂存 basename
              优先, -stub-run 不暂存时回退 manifest 绝对路径。
  --output    输出表路径

DIAMOND 以 `-f 6` 运行 (qseqid sseqid pident length mismatch gapopen
qstart qend sstart send evalue bitscore)。qseqid 即 Prodigal 的 gene id;
本脚本先从 proteins.faa 序列头提取权威 gene id 集合, 再逐行校验 —— 出现
faa 中不存在的 qseqid 直接报错 (键一致性是 Phase 14 总表 join 的前提,
键不一致即验收失败)。每个 gene 只保留最优命中 (bitscore 降序, 并列时
evalue 升序, 再并列时 sseqid 字典序), 保证 join 键唯一。

空 raw (无命中) 与空 faa 均合法: 该 MAG 无行, 输出仅表头。

输出 diamond_hits.tsv:
  meta_id  mag_id  gene  sseqid  pident  length  evalue  bitscore
"""

import argparse
import os
import csv
import sys

OUTPUT_COLUMNS = [
    "meta_id", "mag_id", "gene",
    "sseqid", "pident", "length", "evalue", "bitscore",
]


def resolve(path):
    """解析 manifest 中的文件路径: 任务目录暂存名 (basename) 优先, 回退
    绝对路径。

    真实运行: raw/faa 作为输入暂存进汇总任务目录, basename 即命中。
    -stub-run: 未被 stub 引用的输入不暂存, 回退到 manifest 中的绝对
    路径 (上游任务 work 目录, 文件仍在)。两者都找不到才报错。
    """
    base = os.path.basename(path)
    if os.path.exists(base):
        return base
    if os.path.exists(path):
        return path
    raise SystemExit(
        f"file not found: {path} (also tried staged basename {base})")


def load_manifest(path):
    """读 manifest (meta_id \\t mag_id \\t raw \\t faa, 无表头)。

    返回 [(meta_id, mag_id, raw_path, faa_path)]; 空行跳过 (collectFile
    newLine:true 的尾行), 列数不足 / mag_id 重复直接报错。
    """
    rows = []
    seen_mag = set()
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            fields = line.rstrip("\n").split("\t")
            if not fields or not fields[0]:
                continue
            if len(fields) < 4:
                raise SystemExit(
                    f"manifest {path} line {lineno} malformed "
                    f"(expected 4 columns): {line!r}")
            meta_id, mag_id, raw_path, faa_path = fields[:4]
            if mag_id in seen_mag:
                raise SystemExit(
                    f"manifest {path}: duplicate mag_id {mag_id!r} "
                    f"(line {lineno})")
            seen_mag.add(mag_id)
            rows.append((meta_id, mag_id, resolve(raw_path),
                         resolve(faa_path)))
    return rows


def load_gene_ids(path):
    """从 proteins.faa 序列头提取 gene id 集合 (头行 '>' 后第一个 token)。

    Prodigal 头形如 '>k119_18_1 # 2 # 1482 # 1 # ID=...', gene id 即
    第一个空白分隔 token。首个非空行必须是 FASTA 头 (否则报错), 之后的
    非头行按序列跳过。重复 gene id 报错 (Phase 14 join 键必须唯一)。
    """
    ids = set()
    seen_header = False
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                seen_header = True
                gene = line[1:].split()[0]
                if gene in ids:
                    raise SystemExit(f"{path}: duplicate gene id {gene!r}")
                ids.add(gene)
            elif not seen_header:
                raise SystemExit(
                    f"{path} line {lineno} is not a FASTA header: {line!r}")
    return ids


def _is_better(new, old):
    """new 是否优于 old: bitscore 降序 → evalue 升序 → sseqid 字典序。"""
    new_bits, new_evalue, new_subject = new[0], new[1], new[2]
    old_bits, old_evalue, old_subject = old[0], old[1], old[2]
    if (new_bits, new_evalue) != (old_bits, old_evalue):
        return (new_bits, -new_evalue) > (old_bits, -old_evalue)
    return new_subject < old_subject


def parse_diamond(path, gene_ids):
    """解析一个 DIAMOND outfmt-6 文件, 返回 {gene: (bitscore, evalue,
    sseqid, pident, length)} —— 每 gene 仅最优命中。"""
    best = {}
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                raise SystemExit(
                    f"{path} line {lineno}: expected 12 outfmt-6 columns, "
                    f"got {len(fields)}: {line!r}")
            gene = fields[0]
            if gene not in gene_ids:
                raise SystemExit(
                    f"{path} line {lineno}: unknown gene id {gene!r} "
                    f"(not in proteins.faa)")
            try:
                bitscore = float(fields[11])
                evalue = float(fields[10])
            except ValueError:
                raise SystemExit(
                    f"{path} line {lineno}: non-numeric score columns: "
                    f"{line!r}")
            candidate = (bitscore, evalue, fields[1], fields[2], fields[3])
            prev = best.get(gene)
            if prev is None or _is_better(candidate, prev):
                best[gene] = candidate
    return best


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0])
    parser.add_argument("--manifest", required=True,
                        help="meta_id \\t mag_id \\t raw \\t faa (无表头)")
    parser.add_argument("--output", required=True,
                        help="输出 diamond_hits.tsv 路径")
    args = parser.parse_args(argv)

    rows = []
    for meta_id, mag_id, raw_path, faa_path in load_manifest(args.manifest):
        gene_ids = load_gene_ids(faa_path)
        best = parse_diamond(raw_path, gene_ids)
        for gene in sorted(best):
            bitscore, evalue, sseqid, pident, length = best[gene]
            rows.append([meta_id, mag_id, gene, sseqid, pident, length,
                         evalue, bitscore])
    rows.sort(key=lambda r: (r[0], r[1], r[2]))

    with open(args.output, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_COLUMNS)
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
