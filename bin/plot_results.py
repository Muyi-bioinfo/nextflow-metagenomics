#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""plot_results.py — Phase 21 结果可视化: 读各 Phase 已产出 TSV, 画 7 张图

只消费现成表, 不生产任何矩阵 (矩阵由 Phase 13/20 产出)。技术选型 Python +
matplotlib (Agg, 无显示环境) + numpy —— 不引 R、不引 scipy; PCoA 用 numpy 对
距离矩阵做 Gower 双中心化 + np.linalg.eigh 特征分解。

─── 7 张图 (子命令 → 输出 PNG → 发布目录) ────────────────────────────────────
  mag-abundance-heatmap    mag_abundance_heatmap.png         13_abundance/figures/
  qc-scatter               completeness_vs_contamination.png  08_mag_qc/figures/
  mag-taxonomy-composition mag_taxonomy_composition.png       10_mag_taxonomy/figures/
  taxonomic-composition    taxonomic_composition.png          03_taxonomy/figures/
  beta-pcoa                beta_diversity_pcoa.png            03_taxonomy/figures/
  pathway-heatmap          pathway_abundance_heatmap.png      04_function/figures/
  workflow-summary         workflow_summary.png               99_multiqc/figures/

─── 容忍语义 (与 bin/merge_read_based.py 一致) ────────────────────────────────
  - 空表 (0 字节 / 仅表头 / 0 数据行): 告警跳过, 返回 0 不写 PNG (调用方以
    optional 输出接住) —— stub 与真实 0-MAG / 0-taxa 数据都落在这里。
  - 缺列 (如 mag_qc.tsv 缺 Completeness/Contamination): 明确报错 (SystemExit),
    不静默画错图。
  - 单点 / 少量点: 正常绘制 (散点/漏斗/热图都容忍 n=1)。
  - PCoA: 距离矩阵 <2 样本, 或 Gower 双中心化后正特征值 <2 个 (Bray-Curtis
    非欧氏, 可能只出 1 个正特征值) 时告警跳过。
  - 漏斗: 任一输入层级缺失 (0 字节哨兵 / 0 行) 跳过该级, 画剩余层级; 全缺则
    整体跳过。

依赖: matplotlib + numpy (environment.yml 已追加)。所有子命令写 PNG 用固定
figsize=(8,6)、dpi=100 (输出 800×600, 便于单测断言尺寸)。
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")  # 无显示环境 (headless)
import matplotlib.pyplot as plt
import numpy as np

# 统一图尺寸与 DPI (单测按此断言 PNG 像素尺寸)
FIGSIZE = (8, 6)
DPI = 100


# ---------------------------------------------------------------------------
# 通用工具
# ---------------------------------------------------------------------------

def is_present(path):
    """文件存在且非 0 字节。空表/哨兵 → False。"""
    return bool(path) and os.path.exists(path) and os.path.getsize(path) > 0


def skip(msg):
    """告警跳过 (合法空数据, 非错误)。"""
    print(f"note: {msg}", file=sys.stderr)
    return 0


def save_fig(fig, out_path):
    """固定 dpi 落盘 (无 bbox_inches, 保证尺寸确定)。"""
    fig.tight_layout()
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    print(f"wrote {out_path}", file=sys.stderr)


def load_matrix(path):
    """读宽表矩阵: 首列 = 行标签, 其余列 = 数值样本列。

    返回 (row_labels, col_labels, data: float ndarray); 空表/仅表头 → None。
    缺列行以 0 补齐 (防御性)。"""
    if not is_present(path):
        return None
    with open(path, encoding="utf-8") as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    if len(lines) < 2:
        return None  # 仅表头 (0 数据行)
    header = lines[0].split("\t")
    col_labels = header[1:]
    if not col_labels:
        return None
    row_labels, vals = [], []
    for line in lines[1:]:
        f = line.split("\t")
        row_labels.append(f[0].strip())
        row = f[1:]
        if len(row) < len(col_labels):
            row = row + ["0"] * (len(col_labels) - len(row))
        vals.append([float(x) for x in row[:len(col_labels)]])
    return row_labels, col_labels, np.asarray(vals, dtype=float)


def read_tsv_rows(path):
    """读 TSV → (header, rows); 0 字节 → (None, [])。"""
    if not is_present(path):
        return None, []
    with open(path, encoding="utf-8") as fh:
        first = fh.readline()
        if not first.strip():
            return None, []
        header = first.rstrip("\n").split("\t")
        rows = []
        for line in fh:
            if not line.strip():
                continue
            rows.append(line.rstrip("\n").split("\t"))
    return header, rows


def top_labels(labels, values, top, other="Other"):
    """把 labels/values 按 values 降序取前 top 个, 其余合并为 other。"""
    pairs = sorted(zip(labels, values), key=lambda p: -p[1])
    head = pairs[:top]
    rest = pairs[top:]
    out_labels = [p[0] for p in head]
    out_values = [p[1] for p in head]
    if rest:
        out_labels.append(other)
        out_values.append(sum(p[1] for p in rest))
    return out_labels, out_values


def composition_counts(header, rows, colname):
    """统计某分类列 (如 phylum/class) 各值的非空计数 (MAG 数)。"""
    if header is None:
        return {}
    if colname not in header:
        raise SystemExit(
            f"ERROR: taxonomy table lacks required column {colname!r} "
            f"(got header {header!r})")
    idx = header.index(colname)
    counts = {}
    for row in rows:
        if len(row) <= idx:
            continue
        key = row[idx].strip() or "Unassigned"
        counts[key] = counts.get(key, 0) + 1
    return counts


# ---------------------------------------------------------------------------
# ① 13_abundance — MAG 丰度热图 (本机真实可验证)
# ---------------------------------------------------------------------------

def cmd_abundance_heatmap(args):
    m = load_matrix(args.matrix)
    if m is None:
        return skip("mag_abundance.tsv is empty or has no data rows — skip "
                    "abundance heatmap")
    mags, samples, data = m
    if data.shape[0] == 0 or data.shape[1] == 0:
        return skip("abundance matrix has no MAGs or no samples — skip heatmap")

    fig, ax = plt.subplots(figsize=FIGSIZE)
    im = ax.imshow(data, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(mags)))
    ax.set_yticklabels(mags, fontsize=9)
    # 小矩阵标注数值, 大矩阵省略 (避免拥挤)
    if data.size <= 200:
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                ax.text(j, i, f"{data[i, j]:.2f}", ha="center", va="center",
                        fontsize=7, color="white" if data[i, j] > 0.5 else "black")
    ax.set_xlabel("Sample")
    ax.set_ylabel("MAG")
    ax.set_title("MAG relative abundance")
    fig.colorbar(im, ax=ax, label="Relative abundance (0-1)")
    save_fig(fig, args.output)
    return 0


# ---------------------------------------------------------------------------
# ② 08_mag_qc — 完整度 vs 污染度散点 + 阈值线
# ---------------------------------------------------------------------------

def cmd_qc_scatter(args):
    header, rows = read_tsv_rows(args.qc)
    if header is None or not rows:
        return skip("mag_qc.tsv is empty or has no data rows — skip QC scatter")
    for col in ("Completeness", "Contamination"):
        if col not in header:
            raise SystemExit(
                f"ERROR: mag_qc.tsv lacks required column {col!r} "
                f"(got header {header!r})")
    ci = header.index("Completeness")
    ti = header.index("Contamination")
    comp = [float(r[ci]) for r in rows if len(r) > max(ci, ti)]
    cont = [float(r[ti]) for r in rows if len(r) > max(ci, ti)]

    fig, ax = plt.subplots(figsize=FIGSIZE)
    passed = [c >= args.min_completeness and t <= args.max_contamination
              for c, t in zip(comp, cont)]
    ax.scatter([t for t, p in zip(cont, passed) if p],
               [c for c, p in zip(comp, passed) if p],
               c="tab:green", label="Passed", alpha=0.8, edgecolors="k", s=40)
    ax.scatter([t for t, p in zip(cont, passed) if not p],
               [c for c, p in zip(comp, passed) if not p],
               c="tab:red", label="Failed", alpha=0.8, edgecolors="k", s=40)
    ax.axvline(args.max_contamination, color="tab:red", linestyle="--",
               label=f"Max contamination = {args.max_contamination}")
    ax.axhline(args.min_completeness, color="tab:green", linestyle="--",
               label=f"Min completeness = {args.min_completeness}")
    ax.set_xlabel("Contamination (%)")
    ax.set_ylabel("Completeness (%)")
    ax.set_title("MAG completeness vs contamination")
    ax.legend(fontsize=8)
    save_fig(fig, args.output)
    return 0


# ---------------------------------------------------------------------------
# ③ 10_mag_taxonomy — 门/纲组成 (MAG 计数)
# ---------------------------------------------------------------------------

def cmd_taxonomy_composition(args):
    header, rows = read_tsv_rows(args.taxonomy)
    if header is None or not rows:
        return skip("mag_taxonomy.tsv is empty or has no data rows — skip "
                    "taxonomy composition")
    phylum = composition_counts(header, rows, "phylum")
    klass = composition_counts(header, rows, "class")
    if not phylum and not klass:
        return skip("no phylum/class assignments in mag_taxonomy.tsv — skip")

    fig, axes = plt.subplots(2, 1, figsize=FIGSIZE, sharex=False)
    for ax, counts, title in (
            (axes[0], phylum, "Phylum composition (MAG count)"),
            (axes[1], klass, "Class composition (MAG count)")):
        labels, values = top_labels(list(counts), list(counts.values()), args.top)
        y = np.arange(len(labels))
        ax.barh(y, values, color="tab:blue", alpha=0.85)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9)
        ax.invert_yaxis()  # 最大值在上
        ax.set_xlabel("Number of MAGs")
        ax.set_title(title, fontsize=10)
        for yi, v in zip(y, values):
            ax.text(v, yi, f" {v}", va="center", fontsize=8)
    save_fig(fig, args.output)
    return 0


# ---------------------------------------------------------------------------
# ④ 03_taxonomy — 跨样本 top taxa 组成
# ---------------------------------------------------------------------------

def _pick_level(paths):
    """从 merged_<level>.tsv 文件列表选层级: S 前缀优先, 否则首个。"""
    base = [os.path.basename(p) for p in paths]
    for b in base:
        if b.startswith("merged_") and b[len("merged_"):].startswith("S"):
            return paths[base.index(b)]
    return paths[0]


def cmd_taxonomic_composition(args):
    matrices = list(args.matrices)
    if not matrices:
        return skip("no merged_<level>.tsv matrices provided — skip taxonomic "
                    "composition")
    level_path = _pick_level(matrices)
    m = load_matrix(level_path)
    if m is None:
        return skip(f"{level_path} is empty or has no data rows — skip "
                    "taxonomic composition")
    taxa, samples, data = m
    totals = data.sum(axis=1)  # 跨样本总丰度
    labels, values = top_labels(taxa, list(totals), args.top)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    y = np.arange(len(labels))
    ax.barh(y, values, color="tab:green", alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel(f"Cumulative fraction_total_reads ({len(samples)} samples)")
    ax.set_title(f"Top {min(args.top, len(labels))} taxa (level "
                 f"{os.path.basename(level_path).replace('merged_', '').replace('.tsv', '')})")
    for yi, v in zip(y, values):
        ax.text(v, yi, f" {v:.4g}", va="center", fontsize=8)
    save_fig(fig, args.output)
    return 0


# ---------------------------------------------------------------------------
# ⑤ 03_taxonomy — β 多样性 PCoA (numpy 特征分解)
# ---------------------------------------------------------------------------

def cmd_beta_pcoa(args):
    header, rows = read_tsv_rows(args.distance)
    if header is None or len(rows) < 2:
        return skip("beta_diversity.tsv missing, empty, or <2 samples — skip PCoA")
    samples = header[1:]
    try:
        dist = np.asarray([[float(x) for x in r[1:]] for r in rows], dtype=float)
    except (ValueError, IndexError):
        raise SystemExit(f"ERROR: beta_diversity.tsv malformed distance matrix")
    n = len(samples)
    if dist.shape != (n, n):
        raise SystemExit(
            f"ERROR: beta_diversity.tsv is not a square {n}x{n} matrix "
            f"(got {dist.shape})")

    # Gower 双中心化: B = -0.5 * J * D² * J, J = I - (1/n) 1 1ᵀ
    d2 = dist ** 2
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    pos = eigvals > 1e-10
    if pos.sum() < 2:
        return skip("PCoA yields <2 positive eigenvalues (Bray-Curtis is "
                    "non-Euclidean) — skip")
    coords = eigvecs[:, :2] * np.sqrt(eigvals[:2])
    var = eigvals[:2] / eigvals[pos].sum()

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.scatter(coords[:, 0], coords[:, 1], c="tab:purple", s=60, edgecolors="k")
    for i, s in enumerate(samples):
        ax.annotate(s, (coords[i, 0], coords[i, 1]), fontsize=9,
                    xytext=(4, 4), textcoords="offset points")
    ax.axhline(0, color="grey", linewidth=0.6)
    ax.axvline(0, color="grey", linewidth=0.6)
    ax.set_xlabel(f"PCoA 1 ({100 * var[0]:.1f}% variance)")
    ax.set_ylabel(f"PCoA 2 ({100 * var[1]:.1f}% variance)")
    ax.set_title("Bray-Curtis beta diversity (PCoA)")
    save_fig(fig, args.output)
    return 0


# ---------------------------------------------------------------------------
# ⑥ 04_function — 通路丰度热图 (top-N pathway 按丰度排序)
# ---------------------------------------------------------------------------

def cmd_pathway_heatmap(args):
    m = load_matrix(args.matrix)
    if m is None:
        return skip("merged_pathabundance.tsv is empty or has no data rows — "
                    "skip pathway heatmap")
    pathways, samples, data = m
    totals = data.sum(axis=1)
    order = np.argsort(totals)[::-1][:args.top]
    if order.size == 0:
        return skip("no pathways to plot — skip pathway heatmap")
    sub = data[order]
    sub_paths = [pathways[i] for i in order]
    # 丰度跨度大 (RPK 跨数量级), 用 log10(x+1e-6) 做色阶, 保持 0 → 黑/底端
    display = np.log10(sub + 1e-6)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    im = ax.imshow(display, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(sub_paths)))
    ax.set_yticklabels(sub_paths, fontsize=7)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Pathway")
    ax.set_title(f"Top {len(sub_paths)} pathways by abundance (log10)")
    fig.colorbar(im, ax=ax, label="log10(Abundance)")
    save_fig(fig, args.output)
    return 0


# ---------------------------------------------------------------------------
# ⑦ 99_multiqc — MAG 工作流漏斗
# ---------------------------------------------------------------------------

def _count_raw_bins(path):
    if not is_present(path):
        return False, 0
    n = 0
    with open(path, encoding="utf-8") as fh:
        fh.readline()  # header
        for line in fh:
            if line.strip():
                n += 1
    return n > 0, n


def _count_qc_pass(path, min_c, max_c):
    header, rows = read_tsv_rows(path)
    if header is None or not rows:
        return False, 0
    for col in ("Completeness", "Contamination"):
        if col not in header:
            raise SystemExit(
                f"ERROR: mag_qc.tsv lacks required column {col!r} "
                f"(got header {header!r})")
    ci, ti = header.index("Completeness"), header.index("Contamination")
    n = 0
    for r in rows:
        if len(r) > max(ci, ti) and float(r[ci]) >= min_c and float(r[ti]) <= max_c:
            n += 1
    return True, n


def _count_representatives(path):
    # membership: mag_id \t sample \t representative_mag_id (无表头)
    if not is_present(path):
        return False, 0
    n = 0
    has = False
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            has = True
            if len(f) >= 3 and f[0] == f[2]:
                n += 1
    return has, n


def _count_assigned(path):
    header, rows = read_tsv_rows(path)
    if header is None or not rows:
        return False, 0
    if "domain" not in header:
        raise SystemExit(
            f"ERROR: mag_taxonomy.tsv lacks required column 'domain' "
            f"(got header {header!r})")
    di = header.index("domain")
    n = 0
    for r in rows:
        if len(r) > di and r[di].strip():
            n += 1
    return True, n


def cmd_workflow_summary(args):
    levels = []  # (label, count), 仅保留"有数据"的层级
    raw_ok, raw = _count_raw_bins(args.bin_summary)
    if raw_ok:
        levels.append(("Raw bins", raw))
    qc_ok, qc = _count_qc_pass(args.qc, args.min_completeness, args.max_contamination)
    if qc_ok:
        levels.append(("Passed QC", qc))
    drep_ok, drep = _count_representatives(args.membership)
    if drep_ok:
        levels.append(("After dRep", drep))
    tax_ok, tax = _count_assigned(args.taxonomy)
    if tax_ok:
        levels.append(("GTDB-Tk assigned", tax))
    if not levels:
        return skip("no MAG pipeline tables present — skip workflow summary")

    fig, ax = plt.subplots(figsize=FIGSIZE)
    labels = [l for l, _ in levels]
    counts = [c for _, c in levels]
    y = np.arange(len(labels))
    ax.barh(y, counts, color="tab:orange", alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Number of MAGs")
    ax.set_title("MAG workflow summary")
    for yi, v in zip(y, counts):
        ax.text(v, yi, f" {v}", va="center", fontsize=9)
    save_fig(fig, args.output)
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("mag-abundance-heatmap", help="MAG 丰度热图")
    p.add_argument("--matrix", required=True)
    p.add_argument("--output", required=True)
    p.set_defaults(func=cmd_abundance_heatmap)

    p = sub.add_parser("qc-scatter", help="完整度 vs 污染度散点 + 阈值线")
    p.add_argument("--qc", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--min-completeness", type=float, default=50.0)
    p.add_argument("--max-contamination", type=float, default=10.0)
    p.set_defaults(func=cmd_qc_scatter)

    p = sub.add_parser("mag-taxonomy-composition", help="MAG 门/纲组成")
    p.add_argument("--taxonomy", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--top", type=int, default=10)
    p.set_defaults(func=cmd_taxonomy_composition)

    p = sub.add_parser("taxonomic-composition", help="跨样本 top taxa")
    p.add_argument("--matrices", nargs="+", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--top", type=int, default=20)
    p.set_defaults(func=cmd_taxonomic_composition)

    p = sub.add_parser("beta-pcoa", help="β 多样性 PCoA (numpy 特征分解)")
    p.add_argument("--distance", required=True)
    p.add_argument("--output", required=True)
    p.set_defaults(func=cmd_beta_pcoa)

    p = sub.add_parser("pathway-heatmap", help="top-N 通路丰度热图")
    p.add_argument("--matrix", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--top", type=int, default=50)
    p.set_defaults(func=cmd_pathway_heatmap)

    p = sub.add_parser("workflow-summary", help="MAG 工作流漏斗")
    p.add_argument("--bin-summary", required=True)
    p.add_argument("--qc", required=True)
    p.add_argument("--membership", required=True)
    p.add_argument("--taxonomy", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--min-completeness", type=float, default=50.0)
    p.add_argument("--max-contamination", type=float, default=10.0)
    p.set_defaults(func=cmd_workflow_summary)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
