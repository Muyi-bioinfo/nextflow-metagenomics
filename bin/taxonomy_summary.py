#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Parse GTDB-Tk classify_wf summary files into a per-MAG taxonomy table.

输入:
  --manifest  两列 TSV (mag_id, sample), 由 GTDBTK process 产出 —— 本批次
              提交分类的全部 MAG (去冗余后的代表 MAG) 的权威清单。
  --bac       gtdbtk.bac120.summary.tsv (不存在或为空 = 无细菌基因组)
  --ar        gtdbtk.ar53.summary.tsv   (不存在或为空 = 无古菌基因组)
              两者都不存在/为空也不报错, 只是该域没有分类结果。
  --membership  可选, 三列 TSV (mag_id, sample, representative_mag_id),
              由 dRep 聚类成员表规范化而来, 是**全部 qualified MAG** 的权威
              清单。未提供时 (dRep 未接入) 每个 MAG 视为自身代表, 输出仅由
              manifest + summary 决定。

GTDB-Tk 的 summary user_genome 列即输入 FASTA 的文件名 (<mag_id>.fa)。
剥掉扩展名恢复 mag_id, 再经 manifest 找回 sample —— MAG 身份在批处理
聚合之后依然可追踪。manifest 中的 MAG 若未出现在任何 summary (分类失败),
保留一行空分类, 下游按 mag_id join 时不丢 MAG。

去冗余回填 (--membership): GTDB-Tk 只分类代表 MAG (manifest), 冗余成员
MAG 保留自身身份 (sample / mag_id), 分类继承自同簇代表 MAG, 并在
rep_mag_id 列标明来源; 代表 MAG 该列为空。

输出 mag_taxonomy.tsv:
  sample  mag_id  rep_mag_id  domain  phylum  class  order  family  genus  species

classification 列是 d__;p__;c__;o__;f__;g__;s__ 分号串。GTDB-Tk 在缺级时
保留空 token (如 "...;f__;g__;s__"), 因此按 ';' 以 -1 拆分以免丢尾部空列,
并**按前缀** (d__/p__/...) 归位而非按位置 —— 缺级时不会错位。
"""

import argparse
import csv
import os
import sys

# 常见 FASTA 扩展名, 用于把 user_genome 还原为 MAG ID (长后缀在前)
FASTA_EXTS = (
    ".fa.gz", ".fna.gz", ".fasta.gz", ".faa.gz",
    ".fa", ".fna", ".fasta", ".faa", ".fsa",
)

# classification token 前缀 → 表列名 (缺级时该列为空字符串)
RANK_PREFIXES = (
    ("d__", "domain"),
    ("p__", "phylum"),
    ("c__", "class"),
    ("o__", "order"),
    ("f__", "family"),
    ("g__", "genus"),
    ("s__", "species"),
)
RANK_BY_PREFIX = dict(RANK_PREFIXES)

COLUMNS = ["sample", "mag_id", "rep_mag_id"] + [name for _, name in RANK_PREFIXES]


def load_manifest(path):
    """读 manifest (mag_id \\t sample, 无表头), 返回 {mag_id: sample}。"""
    mapping = {}
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            fields = line.rstrip("\n").split("\t")
            # 空行 (如 collectFile newLine 产生的尾部换行) 直接跳过
            if not fields or not fields[0]:
                continue
            if len(fields) < 2:
                raise SystemExit(f"manifest {path} line {lineno} malformed: {line!r}")
            mag_id, sample = fields[0], fields[1]
            if mag_id in mapping:
                raise SystemExit(f"duplicate mag_id in manifest: {mag_id}")
            mapping[mag_id] = sample
    return mapping


def to_mag_id(user_genome):
    """把 summary 的 user_genome (文件名) 还原为 MAG ID。"""
    for ext in FASTA_EXTS:
        if user_genome.endswith(ext):
            return user_genome[: -len(ext)]
    return user_genome


def load_membership(path):
    """读 dRep 成员表 (mag_id \\t sample \\t representative_mag_id, 无表头)。

    返回 {mag_id: (sample, representative_mag_id)} —— 全部 qualified MAG 的
    权威清单: 代表 MAG 的 representative 指向自身。
    """
    membership = {}
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            fields = line.rstrip("\n").split("\t")
            # 空行 (如 collectFile newLine 产生的尾部换行) 直接跳过
            if not fields or not fields[0]:
                continue
            if len(fields) < 3:
                raise SystemExit(f"membership {path} line {lineno} malformed: {line!r}")
            mag_id, sample, rep = fields[0], fields[1], fields[2]
            if mag_id in membership:
                raise SystemExit(f"duplicate mag_id in membership: {mag_id}")
            membership[mag_id] = (sample, rep)
    return membership


def parse_classification(raw, user_genome):
    """拆 d__;p__;c__;o__;f__;g__;s__ 串; 缺级 token 为空, 按前缀归位。"""
    levels = {name: "" for _, name in RANK_PREFIXES}
    classification = raw.strip()
    if classification.startswith("Unclassified"):
        return levels
    for token in classification.split(";"):
        token = token.strip()
        if not token:
            continue
        rank = RANK_BY_PREFIX.get(token[:3])
        if rank is None:
            raise SystemExit(
                f"unrecognised rank token {token!r} in classification of "
                f"{user_genome}: {classification!r}"
            )
        levels[rank] = token[3:]
    return levels


def read_summary(path, manifest, rows_by_mag, seen):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        print(f"note: {path or '(none)'} absent or empty — skipped", file=sys.stderr)
        return
    with open(path, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for required in ("user_genome", "classification"):
            if required not in reader.fieldnames:
                raise SystemExit(f"{path} lacks required column {required!r}")
        for row in reader:
            user_genome = row["user_genome"].strip()
            mag_id = to_mag_id(user_genome)
            if mag_id not in manifest:
                raise SystemExit(
                    f"user_genome {user_genome!r} (→ mag_id {mag_id!r}) not found "
                    f"in manifest — MAG ID mapping broken"
                )
            if mag_id in seen:
                raise SystemExit(f"duplicate user_genome mapping to mag_id {mag_id!r}")
            seen.add(mag_id)
            levels = parse_classification(row["classification"], user_genome)
            rows_by_mag[mag_id] = [manifest[mag_id], mag_id] + [
                levels[name] for _, name in RANK_PREFIXES
            ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--membership", default="")
    parser.add_argument("--bac", default="")
    parser.add_argument("--ar", default="")
    args = parser.parse_args()

    manifest = load_manifest(args.manifest)
    membership = load_membership(args.membership) if args.membership else None
    if membership is not None:
        # 成员表是全部 qualified MAG 的权威清单, manifest (GTDB-Tk 输入,
        # 通常为去冗余后的代表 MAG) 必须完整落在其中
        for mag_id in manifest:
            if mag_id not in membership:
                raise SystemExit(
                    f"manifest mag {mag_id!r} missing from membership — wiring broken"
                )

    rows_by_mag, seen = {}, set()
    read_summary(args.bac, manifest, rows_by_mag, seen)
    read_summary(args.ar, manifest, rows_by_mag, seen)

    # manifest 中的 MAG 若未出现在任何 summary (分类失败), 保留一行空分类
    for mag_id, sample in manifest.items():
        if mag_id not in rows_by_mag:
            rows_by_mag[mag_id] = [sample, mag_id] + [""] * len(RANK_PREFIXES)

    if membership is not None:
        # 去冗余回填: 代表 MAG 用自身分类 (rep_mag_id 留空); 冗余成员 MAG
        # 保留自身身份 (sample / mag_id), 分类继承自同簇代表 MAG。
        for mag_id, (sample, rep) in membership.items():
            if rep == mag_id:
                if mag_id not in manifest:
                    raise SystemExit(
                        f"representative {mag_id!r} not in GTDB-Tk manifest — wiring broken"
                    )
                row = rows_by_mag[mag_id]
                rows_by_mag[mag_id] = row[:2] + [""] + row[2:]
            else:
                if rep not in manifest:
                    raise SystemExit(
                        f"representative {rep!r} of {mag_id!r} not in "
                        f"GTDB-Tk manifest — mapping broken"
                    )
                if mag_id in manifest:
                    raise SystemExit(
                        f"{mag_id!r} is both in GTDB-Tk manifest and a "
                        f"non-representative in membership — inconsistent input"
                    )
                rep_row = rows_by_mag[rep]
                rows_by_mag[mag_id] = [sample, mag_id, rep] + rep_row[3:]
    else:
        # 无成员表 (dRep 未接入): 每个 MAG 视为自身代表, rep_mag_id 列留空
        for mag_id in rows_by_mag:
            row = rows_by_mag[mag_id]
            rows_by_mag[mag_id] = row[:2] + [""] + row[2:]

    rows = sorted(rows_by_mag.values(), key=lambda r: (r[0], r[1]))
    with open(args.output, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        writer.writerows(rows)


if __name__ == "__main__":
    main()
