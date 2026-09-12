#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""merge_read_based.py — Phase 20 read-based 跨样本整合: 长表 → 宽表矩阵

把逐样本的 read-based 结果合并成跨样本矩阵, 为跨样本比较 (PCoA / 丰度热图 /
样本聚类) 提供数据源。本脚本**只产出数据矩阵, 不画图**。

两个子命令:
  bracken        Bracken 各层级丰度表 → 样本×taxa 矩阵 (每个层级一张)
                 + (可选) Bray-Curtis 距离矩阵 (样本×样本)
  pathabundance  HUMAnN pathabundance → 样本×pathway 矩阵

─── 数据口径 ────────────────────────────────────────────────────────────────
  Bracken 合并值 = 逐样本表的 `fraction_total_reads` 列 (0-1 相对丰度, 已按
    Bracken 输出的原始字符串原样保留, 不做重归一化); 行键 = `name` 列
    (分类单元名)。某样本未检出的分类单元补 "0"。
  HUMAnN pathabundance 合并值 = 逐样本表的 Abundance 列 (RPK 通路丰度);
    行键 = pathway 字符串 (含 "id: name" 前缀, 原样)。缺失 pathway/样本补 "0"。
  Bray-Curtis 基于所选层级 (默认 S 前缀层级, 否则第一个层级) 的相对丰度向量;
    距离 = Σ|a-b| / Σ(a+b), 对角线 0, 对称。仅纯 Python 标准库实现 —— 不引入
    numpy/scipy (与 bin/ 其余脚本一致, 便于 conda/容器/单测三处复用)。

─── 路径解析 (resolve 双路径) ───────────────────────────────────────────────
  manifest 中记录的是上游 work 目录的绝对路径 (collectFile `${path}` 插值,
  Phase 12/14 模式)。真实运行时任务目录暂存 basename 命中优先; -stub-run 未
  被 stub 引用的输入不暂存, 回退 manifest 绝对路径 (上游 work 文件仍在)。
  两者都找不到才报错。

─── 容忍语义 ────────────────────────────────────────────────────────────────
  - 空表 (0 字节 / 仅表头): 该样本贡献 0 行, 不报错 —— stub 与真实 0-read
    样本都落到这里 (不伪造)。
  - 缺列: 缺少 `name` / `fraction_total_reads` (bracken) 或 pathway 两列
    (pathabundance) 时明确报错 (接线/格式错误)。
  - 单样本 / 0 taxa: Bray-Curtis 无意义, 告警并跳过 (不产出距离矩阵, 调用方
    以 optional 输出接住)。
  - 同一分类单元名在同一文件内重复: 报错 (键冲突, 不静默覆盖)。

仅依赖标准库。
"""
import argparse
import csv
import os
import sys


def resolve(path):
    """manifest 中的文件路径 → 可打开路径: 任务目录 basename 优先, 回退
    绝对路径 (Phase 12/14 双路径模式)。"""
    base = os.path.basename(path)
    if os.path.exists(base):
        return base
    if os.path.exists(path):
        return path
    raise SystemExit(f"ERROR: cannot resolve {path!r}: neither staged basename "
                     f"nor manifest absolute path exists")


def read_manifest(path, ncols, what):
    """读 manifest (无表头 TSV)。bracken: 3 列 (level/sample/file);
    pathabundance: 2 列 (sample/file)。返回行列表 (按文件顺序, manifest 已
    collectFile sort:true 排序)。"""
    rows = []
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < ncols:
                raise SystemExit(f"manifest {path} line {lineno} malformed "
                                 f"(expect {ncols} cols): {line!r}")
            rows.append(tuple(fields[:ncols]))
    if not rows:
        raise SystemExit(f"ERROR: {what} manifest {path} is empty — no per-sample "
                         f"tables to merge")
    return rows


# ---------------------------------------------------------------------------
# 逐样本表解析
# ---------------------------------------------------------------------------

def parse_bracken(path):
    """Bracken 丰度表 → [(name, fraction_total_reads), ...] (文件顺序)。
    空表 (0 字节 / 无表头) → []; 缺列报错。"""
    p = resolve(path)
    if os.path.getsize(p) == 0:
        return []
    with open(p, encoding="utf-8") as fh:
        first = fh.readline()
        if not first.strip():
            return []
        header = first.rstrip("\n").split("\t")
        for col in ("name", "fraction_total_reads"):
            if col not in header:
                raise SystemExit(f"ERROR: bracken table {path} lacks required "
                                 f"column {col!r} (got header {header!r})")
        name_i = header.index("name")
        frac_i = header.index("fraction_total_reads")
        rows = []
        for lineno, line in enumerate(fh, 2):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(name_i, frac_i):
                raise SystemExit(f"ERROR: bracken table {path} line {lineno} "
                                 f"has too few columns: {line!r}")
            name = fields[name_i].strip()
            frac = fields[frac_i].strip()
            if not name:
                raise SystemExit(f"ERROR: bracken table {path} line {lineno} "
                                 f"has empty name")
            rows.append((name, frac))
    return rows


def parse_pathabundance(path):
    """HUMAnN pathabundance → [(pathway, abundance), ...] (文件顺序)。
    跳过 '#' 注释头; 空表 → []; 缺列报错。"""
    p = resolve(path)
    if os.path.getsize(p) == 0:
        return []
    rows = []
    with open(p, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip():
                continue
            if line.startswith("#"):
                continue  # 注释头 "# Pathway\tAbundance"
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise SystemExit(f"ERROR: pathabundance table {path} line "
                                 f"{lineno} has too few columns: {line!r}")
            pathway = fields[0].strip()
            abundance = fields[1].strip()
            if not pathway:
                raise SystemExit(f"ERROR: pathabundance table {path} line "
                                 f"{lineno} has empty pathway")
            rows.append((pathway, abundance))
    return rows


# ---------------------------------------------------------------------------
# 宽表构建
# ---------------------------------------------------------------------------

def _accumulate(manifest_rows, parser):
    """manifest 行 → {sample: {key: value}}, 保留样本/键的有序并集。"""
    samples = []
    seen_samples = set()
    data = {}
    for row in manifest_rows:
        sample = row[0]
        path = row[-1]
        if sample not in seen_samples:
            seen_samples.add(sample)
            samples.append(sample)
        values = {}
        for key, val in parser(path):
            if key in values:
                raise SystemExit(f"ERROR: duplicate key {key!r} in {path} — "
                                 f"ambiguous, cannot merge")
            values[key] = val
        data[sample] = values
    return samples, data


def write_wide_matrix(out_path, row_header, rows, samples, data):
    """写样本×特征宽表: 首列 row_header, 其余列 = samples (有序), 缺失补 0。"""
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow([row_header] + samples)
        for key in rows:
            writer.writerow([key] + [data[s].get(key, "0") for s in samples])


def union_keys(data):
    """全部样本键的并集, 排序保证输出确定 (-resume 哈希稳定)。"""
    keys = set()
    for values in data.values():
        keys.update(values)
    return sorted(keys)


# ---------------------------------------------------------------------------
# Bray-Curtis 距离
# ---------------------------------------------------------------------------

def bray_curtis(samples, data, taxa):
    """基于相对丰度向量算 Bray-Curtis 距离矩阵 (对称, 对角线 0)。
    返回 [[float]] 形状 len(samples) × len(samples)。"""
    n = len(samples)
    vecs = [[float(data[s].get(t, "0")) for t in taxa] for s in samples]
    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            a, b = vecs[i], vecs[j]
            num = sum(abs(x - y) for x, y in zip(a, b))
            den = sum(x + y for x, y in zip(a, b))
            d = (num / den) if den > 0 else 0.0
            dist[i][j] = dist[j][i] = d
    return dist


def write_distance_matrix(out_path, samples, dist):
    """写样本×样本距离矩阵: 首行首格为空 (标准距离矩阵格式), 行标签 = 样本。"""
    n = len(samples)
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow([""] + samples)
        for i, s in enumerate(samples):
            writer.writerow([s] + [f"{dist[i][j]:.8g}" for j in range(n)])


# ---------------------------------------------------------------------------
# 子命令
# ---------------------------------------------------------------------------

def _pick_beta_level(levels):
    """Beta diversity 层级: 首选 S 前缀层级 (物种族), 否则第一个层级。"""
    for lvl in levels:
        if lvl.startswith("S"):
            return lvl
    return levels[0]


def cmd_bracken(args):
    manifest = read_manifest(args.manifest, 3, "bracken")
    # 层级有序并集 (manifest 已排序, 去重保序)
    levels = []
    seen = set()
    for level, _sample, _path in manifest:
        if level not in seen:
            seen.add(level)
            levels.append(level)

    for level in levels:
        # 该层级的 (sample, path) 行, 样本有序
        rows = [(s, p) for lv, s, p in manifest if lv == level]
        samples, data = _accumulate(rows, parse_bracken)
        keys = union_keys(data)
        out = os.path.join(args.output_dir, f"merged_{level}.tsv")
        write_wide_matrix(out, "name", keys, samples, data)
        print(f"wrote {out} ({len(keys)} taxa × {len(samples)} samples)",
              file=sys.stderr)

    # Beta diversity: 相对丰度向量 (用所选层级的矩阵口径, 缺失补 0)
    beta_level = _pick_beta_level(levels)
    rows = [(s, p) for lv, s, p in manifest if lv == beta_level]
    samples, data = _accumulate(rows, parse_bracken)
    keys = union_keys(data)
    if len(samples) < 2 or not keys:
        print(f"note: skipping Bray-Curtis distance — {len(samples)} sample(s) / "
              f"{len(keys)} taxa at level {beta_level!r} (need ≥2 samples and ≥1 "
              f"taxon)", file=sys.stderr)
        return 0
    dist = bray_curtis(samples, data, keys)
    out = os.path.join(args.output_dir, "beta_diversity.tsv")
    write_distance_matrix(out, samples, dist)
    print(f"wrote {out} ({len(samples)} × {len(samples)}, level {beta_level!r})",
          file=sys.stderr)
    return 0


def cmd_pathabundance(args):
    manifest = read_manifest(args.manifest, 2, "pathabundance")
    samples, data = _accumulate(
        [(s, p) for s, p in manifest], parse_pathabundance)
    keys = union_keys(data)
    write_wide_matrix(args.output, "pathway", keys, samples, data)
    print(f"wrote {args.output} ({len(keys)} pathways × {len(samples)} samples)",
          file=sys.stderr)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_bracken = sub.add_parser("bracken", help="Bracken 丰度 → 样本×taxa 矩阵")
    p_bracken.add_argument("--manifest", required=True,
                           help="3 列 TSV: level \\t sample \\t file (无表头)")
    p_bracken.add_argument("--output-dir", default=".",
                           help="写 merged_<level>.tsv 与 beta_diversity.tsv 的目录")
    p_bracken.set_defaults(func=cmd_bracken)

    p_humann = sub.add_parser("pathabundance", help="HUMAnN pathabundance → 样本×pathway 矩阵")
    p_humann.add_argument("--manifest", required=True,
                          help="2 列 TSV: sample \\t file (无表头)")
    p_humann.add_argument("--output", default="merged_pathabundance.tsv")
    p_humann.set_defaults(func=cmd_pathabundance)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
