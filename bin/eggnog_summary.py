#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Parse per-MAG eggNOG-mapper annotations into one keyed COG/GO/EC/KO table.

输入:
  --manifest  四列 TSV (meta_id, mag_id, emapper_annotations, proteins_faa),
              无表头, 由 ANNOTATION 子工作流 collectFile 产出。文件路径经 resolve()
              解析: 真实运行的暂存 basename 优先, -stub-run 不暂存时回退
              manifest 绝对路径。
  --output    输出表路径

emapper.py 2.1.x 的 .emapper.annotations 首行以 '#' 开头且为列名行
(#query seed_ortholog ... COG_category Description Preferred_name GOs EC
KEGG_ko KEGG_Pathway ...)。列序随版本可能变动, 因此**按列名定位**而非按
位置: #query 列缺失直接报错, 其余注释列缺失时整列留空 (降级不报错)。
查询名与 proteins.faa 的 gene id 逐行校验, 未知 gene id 报错; 同一
query 出现多行报错 (Phase 14 join 键必须唯一)。

eggNOG 的 '-' (无注释占位) 规范化为空字符串; KO / GO 的多值逗号串与
COG 类别字母**原样保留**, 拆分交给 Phase 14。空 annotations (仅表头)
合法: 该 MAG 无行, 输出仅表头。

输出 eggnog_annotations.tsv:
  meta_id  mag_id  gene  COG_category  Description  Preferred_name  EC  KO  GO  KEGG_Pathway
"""

import argparse
import os
import csv
import sys

# emapper annotations 表头列名 → 输出列名 ('#query' 必选, 其余可选)
ANNOTATION_COLUMNS = (
    ("#query", "gene"),
    ("COG_category", "COG_category"),
    ("Description", "Description"),
    ("Preferred_name", "Preferred_name"),
    ("EC", "EC"),
    ("KEGG_ko", "KO"),
    ("GOs", "GO"),
    ("KEGG_Pathway", "KEGG_Pathway"),
)

OUTPUT_COLUMNS = ["meta_id", "mag_id", "gene", "COG_category",
                  "Description", "Preferred_name", "EC", "KO", "GO",
                  "KEGG_Pathway"]


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


def _norm(value):
    """eggNOG 缺失值规范化: '-' (无注释占位) → 空字符串, 其余原样保留
    (KO/GO 逗号串、COG 类别字母不做拆分)。"""
    if value is None:
        return ""
    value = str(value)
    return "" if value == "-" else value


def parse_annotations(path, gene_ids):
    """解析一个 .emapper.annotations 文件, 返回 {gene: [注释字段...]}。

    返回的列表按 OUTPUT_COLUMNS 去掉 meta_id/mag_id/gene 后的顺序排列。
    空文件 (0 字节, 理论上的边界情况) 视为无注释。
    """
    by_name = {}
    rows = {}
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if lineno == 1 or by_name == {}:
                if not line.startswith("#"):
                    raise SystemExit(
                        f"{path} line {lineno}: expected '#'-prefixed "
                        f"header row, got: {line!r}")
                # 列名同时以原始形式与剥 '#' 形式登记, 兼容 '#' 前缀缺失
                raw_header = line.split("\t")
                by_name = {name: idx for idx, name in enumerate(raw_header)}
                by_name.update({name.lstrip("#"): idx
                                for idx, name in enumerate(raw_header)})
                if "#query" not in by_name:
                    raise SystemExit(
                        f"{path}: header missing required '#query' column")
                continue
            fields = line.split("\t")

            def get(name):
                idx = by_name.get(name)
                if idx is None or idx >= len(fields):
                    return ""
                return _norm(fields[idx])

            gene = get("#query")
            if not gene:
                raise SystemExit(
                    f"{path} line {lineno}: empty '#query' value")
            if gene not in gene_ids:
                raise SystemExit(
                    f"{path} line {lineno}: unknown gene id {gene!r} "
                    f"(not in proteins.faa)")
            if gene in rows:
                raise SystemExit(
                    f"{path} line {lineno}: duplicate gene id {gene!r}")
            rows[gene] = [get(name) for name, _ in ANNOTATION_COLUMNS[1:]]
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0])
    parser.add_argument("--manifest", required=True,
                        help="meta_id \\t mag_id \\t raw \\t faa (无表头)")
    parser.add_argument("--output", required=True,
                        help="输出 eggnog_annotations.tsv 路径")
    args = parser.parse_args(argv)

    rows = []
    for meta_id, mag_id, raw_path, faa_path in load_manifest(args.manifest):
        gene_ids = load_gene_ids(faa_path)
        annotations = parse_annotations(raw_path, gene_ids)
        for gene in sorted(annotations):
            rows.append([meta_id, mag_id, gene] + annotations[gene])
    rows.sort(key=lambda r: (r[0], r[1], r[2]))

    with open(args.output, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_COLUMNS)
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
