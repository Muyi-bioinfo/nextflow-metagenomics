#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Parse per-MAG RGI (CARD) JSON outputs into one keyed ARG table.

输入:
  --manifest  四列 TSV (meta_id, mag_id, rgi_json, proteins_faa), 无表头,
              由 ANNOTATION 子工作流 collectFile 产出。文件路径经 resolve()
              解析: 真实运行的暂存 basename 优先, -stub-run 不暂存时回退
              manifest 绝对路径。
  --output    输出表路径

`rgi main -t protein` 的输出 <prefix>.json: 顶层 dict, 键为序列 id, 值为
该命中的注释 dict (含 ORF_ID / ARO_name / AMR_Gene_Family 等字段)。蛋白
模式下 ORF_ID 即输入 proteins.faa 的序列头 = Prodigal gene id; 仍逐行与
faa 头集合校验, 未知 gene id 报错 (键一致性是 Phase 14 join 的前提)。
兼容两种结构: dict 键控 (缺 ORF_ID 字段时回退用 dict 键) 与 list。

同 gene 多条命中 (不同模型/阈值) 按 pct_identity 降序、ARO_name 字典序
去重为一行, 保证 Phase 14 join 键唯一。空文件 (0 字节) 与空 dict {} 均
视为无命中: 该 MAG 无行, 输出仅表头。非空但解析失败 (非法 JSON / 顶层
非 dict 非 list) 报错。

输出 rgi_annotations.tsv:
  meta_id  mag_id  gene  ARO  ARO_accession  AMR_gene_family  drug_class
  resistance_mechanism  pct_identity  model_type
"""

import argparse
import os
import csv
import json
import sys

# RGI JSON 字段名 → 输出列名 (字段缺失留空, 不报错)
JSON_FIELDS = (
    ("ARO_name", "ARO"),
    ("ARO_accession", "ARO_accession"),
    ("AMR_Gene_Family", "AMR_gene_family"),
    ("Drug_Class", "drug_class"),
    ("Resistance_Mechanism", "resistance_mechanism"),
    ("Model_type", "model_type"),
)

OUTPUT_COLUMNS = ["meta_id", "mag_id", "gene", "ARO", "ARO_accession",
                  "AMR_gene_family", "drug_class",
                  "resistance_mechanism", "pct_identity", "model_type"]


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
    """null / None → 空字符串, 其余转字符串原样保留。"""
    if value is None:
        return ""
    return "" if str(value) in ("None", "null") else str(value)


def _to_float(value):
    """百分比转 float (仅用于排序), 无法解析时按 0 处理。"""
    try:
        return float(str(value).replace("%", ""))
    except (TypeError, ValueError):
        return 0.0


def _pct_raw(entry):
    """原始百分比字符串: Percentage_Identity 优先, 回退 Best_Identities。"""
    for key in ("Percentage_Identity", "Best_Identities"):
        if entry.get(key) is not None:
            return _norm(entry.get(key))
    return ""


def _is_better(new, old):
    """new 是否优于 old: pct_identity 降序 → ARO_name 字典序。"""
    new_pct, new_aro = new
    old_pct, old_aro = old
    if new_pct != old_pct:
        return new_pct > old_pct
    return new_aro < old_aro


def parse_rgi(path, gene_ids):
    """解析一个 RGI JSON 文件, 返回 {gene: {field: value, ...}} ——
    每 gene 仅保留最优命中条目 (含 'pct_identity' 与 'ARO' 原始值)。"""
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    if not text.strip():
        return {}
    try:
        data = json.loads(text)
    except json.JSONDecodeError as exc:
        raise SystemExit(f"{path}: invalid RGI JSON: {exc}")

    entries = []
    if isinstance(data, dict):
        for key, value in data.items():
            if isinstance(value, dict):
                entry = dict(value)
                entry.setdefault("ORF_ID", str(key))
                entries.append(entry)
    elif isinstance(data, list):
        for value in data:
            if isinstance(value, dict):
                entries.append(dict(value))
    else:
        raise SystemExit(
            f"{path}: unsupported RGI JSON structure "
            f"(top-level must be dict or list, got {type(data).__name__})")

    best = {}
    for entry in entries:
        gene = _norm(entry.get("ORF_ID"))
        if not gene:
            raise SystemExit(
                f"{path}: entry without ORF_ID / sequence key: {entry!r}")
        if gene not in gene_ids:
            raise SystemExit(
                f"{path}: unknown gene id {gene!r} "
                f"(not in proteins.faa)")
        pct_raw = _pct_raw(entry)
        candidate = (_to_float(pct_raw), _norm(entry.get("ARO_name")))
        prev = best.get(gene)
        if prev is None or _is_better(candidate, prev[0]):
            best[gene] = (candidate, entry, pct_raw)
    return best


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0])
    parser.add_argument("--manifest", required=True,
                        help="meta_id \\t mag_id \\t raw \\t faa (无表头)")
    parser.add_argument("--output", required=True,
                        help="输出 rgi_annotations.tsv 路径")
    args = parser.parse_args(argv)

    rows = []
    for meta_id, mag_id, raw_path, faa_path in load_manifest(args.manifest):
        gene_ids = load_gene_ids(faa_path)
        best = parse_rgi(raw_path, gene_ids)
        for gene in sorted(best):
            _, entry, pct_raw = best[gene]
            aro, aro_acc, family, drug, mech, model = [
                _norm(entry.get(name)) for name, _ in JSON_FIELDS]
            rows.append([meta_id, mag_id, gene, aro, aro_acc, family,
                         drug, mech, pct_raw, model])
    rows.sort(key=lambda r: (r[0], r[1], r[2]))

    with open(args.output, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_COLUMNS)
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
