#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""integrate_functional.py — Phase 14 整合: 三张注释表 → mag_functional_annotation.tsv

输入:
  --manifest  三列 TSV (meta_id, mag_id, fasta_path), 无表头 —— 代表 MAG 的
              权威清单: 用于校验注释表行的 MAG 身份 (未知 mag_id 直接报错,
              键一致性是总表 join 的前提), 并在三表全缺时作为调度触发
              (输出仅表头)。
  --diamond   可选: diamond_hits.tsv (meta_id/mag_id/gene/sseqid/...)。
              本表只提供 (meta_id, mag_id, gene) 键并集, 命中细节
              (sseqid/pident/evalue/bitscore) 不进目标列。
  --eggnog    可选: eggnog_annotations.tsv (meta_id/mag_id/gene/COG_category/
              Description/.../KO/GO/KEGG_Pathway)。
  --rgi       可选: rgi_annotations.tsv (meta_id/mag_id/gene/ARO/...)。
  --pathway   可选: KO→pathway 两列映射 (无表头), 即 params.pathway_db;
              缺失或 0 字节时 Pathway 列留空 (不伪造)。KO 列容忍 "ko:" 前缀
              (加载时剥离), pathway 原样保留; 同一 KO 可多行映射多个 pathway。

三张表按 (meta_id, mag_id, gene) 外连接: 行 = 三表键的并集, 缺失表/缺失
字段留空; 同一键出现在多张表时合并为一行 (同键 meta_id 不一致报错)。任一
表 0 字节视为"该注释分支未运行", 与"表存在但无行"区分 (后者合法, 贡献零行)。

目标列 (固定列序):
  MAG_ID / Gene / KO / COG / GO / Pathway / ARG
    KO      ← eggnog KO (多值逗号串原样保留)
    COG     ← eggnog COG_category (类别字母串原样保留, 如 "EG")
    GO      ← eggnog GO (多值逗号串原样保留)
    Pathway ← KO 各 token 在 pathway_db 中的映射并集, 去重 (保持 KO 内
              token 顺序) 后分号连接; pathway_db 缺失时留空
    ARG     ← rgi ARO (命中名称; ARO_accession/drug_class 等细节不进目标列)

输出按 (MAG_ID, Gene) 排序 (同键再按 meta_id), 内容确定。
"""
import argparse
import csv
import os
import sys

OUTPUT_COLUMNS = ["MAG_ID", "Gene", "KO", "COG", "GO", "Pathway", "ARG"]


def read_table(path, column_pairs):
    """读一张键控注释表 → {key: {目标列: 值}}; 0 字节视为表不存在 (None)。

    column_pairs = [(源列名, 目标列名), ...]; 键列 (meta_id/mag_id/gene)
    固定必选。"""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        print(f"note: annotation table {path or '(none)'} absent or empty — skipped",
              file=sys.stderr)
        return None
    table = {}
    with open(path, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for required in ["meta_id", "mag_id", "gene"] + [src for src, _ in column_pairs]:
            if required not in reader.fieldnames:
                raise SystemExit(f"annotation table {path} lacks required "
                                 f"column {required!r}")
        for row in reader:
            key = (row["meta_id"].strip(), row["mag_id"].strip(), row["gene"].strip())
            if key in table:
                raise SystemExit(f"duplicate key {key!r} in annotation table {path}")
            table[key] = {target: row.get(src, "").strip() for src, target in column_pairs}
    return table


def load_pathway_map(path):
    """读 KO→pathway 映射 (两列, 无表头) → {KO: [pathway, ...]}。

    KO 列剥离 "ko:" 前缀 (KEGG ko_pathway 惯例); 映射行顺序保留, 同一
    KO 多行 = 多个 pathway。0 字节视为未提供 (Pathway 列留空)。
    """
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        print(f"note: pathway mapping {path or '(none)'} absent or empty — "
              f"Pathway column blank", file=sys.stderr)
        return {}
    mapping = {}
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise SystemExit(f"pathway mapping {path} line {lineno} malformed: {line!r}")
            ko = fields[0].strip()
            if ko.startswith("ko:"):
                ko = ko[3:]
            pathway = fields[1].strip()
            if not ko or not pathway:
                raise SystemExit(f"pathway mapping {path} line {lineno}: empty "
                                 f"KO or pathway: {line!r}")
            mapping.setdefault(ko, []).append(pathway)
    return mapping


def ko_to_pathways(ko_string, mapping):
    """KO 多值逗号串 → pathway 并集 (保持 KO 内 token 顺序去重, 分号连接)。"""
    seen, pathways = set(), []
    for token in ko_string.split(","):
        token = token.strip()
        if not token:
            continue
        for pathway in mapping.get(token, []):
            if pathway not in seen:
                seen.add(pathway)
                pathways.append(pathway)
    return ";".join(pathways)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--diamond", default="")
    ap.add_argument("--eggnog", default="")
    ap.add_argument("--rgi", default="")
    ap.add_argument("--pathway", default="")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    # 代表 MAG 权威清单 (fasta_path 列本脚本不用, 仅作调度触发与身份校验)
    manifest = {}
    with open(args.manifest, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            fields = line.rstrip("\n").split("\t")
            if not fields or not fields[0]:
                continue
            if len(fields) < 3:
                raise SystemExit(f"manifest {args.manifest} line {lineno} malformed: {line!r}")
            mag_id = fields[1]
            if mag_id in manifest:
                raise SystemExit(f"duplicate mag_id in manifest: {mag_id}")
            manifest[mag_id] = fields[0]

    diamond = read_table(args.diamond, [])
    eggnog = read_table(args.eggnog, [("COG_category", "COG"), ("KO", "KO"),
                                      ("GO", "GO")])
    rgi = read_table(args.rgi, [("ARO", "ARO")])
    pathway_map = load_pathway_map(args.pathway)

    # (meta_id, mag_id, gene) 外连接: 键并集, 跨表合并为一行
    merged = {}
    for table in (diamond, eggnog, rgi):
        if table is None:
            continue
        for key, fields in table.items():
            meta_id, mag_id, gene = key
            if mag_id not in manifest:
                raise SystemExit(f"annotation row {key!r} has mag_id not in "
                                 f"representative manifest — wiring broken")
            if key not in merged:
                merged[key] = {"meta_id": meta_id, "mag_id": mag_id, "gene": gene,
                               "KO": "", "COG": "", "GO": "", "ARO": ""}
            row = merged[key]
            for column, value in fields.items():
                if row[column] and value and row[column] != value:
                    raise SystemExit(f"key {key!r} has conflicting values for "
                                     f"column {column!r}")
                if value:
                    row[column] = value

    # 键一致性: 同一 (mag_id, gene) 不得出现于两个 meta_id 之下 (Phase 12
    # 契约: 基因 id 按 MAG 生成, 跨 meta 重名即接线断裂)
    seen_pairs = {}
    for meta_id, mag_id, gene in merged:
        pair = (mag_id, gene)
        if pair in seen_pairs and seen_pairs[pair] != meta_id:
            raise SystemExit(f"(mag_id, gene) {pair!r} appears under two meta_ids: "
                             f"{seen_pairs[pair]!r} and {meta_id!r} — wiring broken")
        seen_pairs[pair] = meta_id

    rows = []
    for key in sorted(merged, key=lambda k: (k[1], k[2], k[0])):
        row = merged[key]
        rows.append([
            row["mag_id"],
            row["gene"],
            row["KO"],
            row["COG"],
            row["GO"],
            ko_to_pathways(row["KO"], pathway_map),
            row["ARO"],
        ])

    with open(args.output, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_COLUMNS)
        writer.writerows(rows)


if __name__ == "__main__":
    main()
