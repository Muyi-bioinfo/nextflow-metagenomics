#!/usr/bin/env python3
"""Render pipeline workflow diagrams as PNGs.

Usage:
    conda activate nf-meta
    python3 bin/draw_workflow.py --layout vertical   # docs/images/workflow.png
    python3 bin/draw_workflow.py --layout compact    # docs/images/workflow_compact.png

Drawn with matplotlib (see environment.yml); Agg backend, no display required.
Boxes are sized to their text so long tool chains never overflow.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

REPO_ROOT = Path(__file__).resolve().parent.parent

COLORS = {
    "process": dict(fc="#ffffff", ec="#34495e", tc="#1c2833"),
    "data":    dict(fc="#f2f4f4", ec="#7f8c8d", tc="#1c2833"),
    "output":  dict(fc="#eaf2fb", ec="#2980b9", tc="#1a5276"),
}
ARROW = "#5d6d7e"


# ---------------------------------------------------------------------------
# shared drawing helpers (work in whatever coordinate system the boxes use)
# ---------------------------------------------------------------------------
def anchor(box, side, invert_y=False):
    x, y, w, h = box["x"], box["y"], box["w"], box["h"]
    # vertical layout uses a normal (upward) y-axis; compact uses an inverted
    # (downward) y-axis, which flips which of y±h/2 is the visual top.
    top = (y - h / 2) if invert_y else (y + h / 2)
    bottom = (y + h / 2) if invert_y else (y - h / 2)
    return {
        "top":    (x, top),
        "bottom": (x, bottom),
        "left":   (x - w / 2, y),
        "right":  (x + w / 2, y),
    }[side]


def draw_arrow(ax, pts):
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    ax.plot(xs, ys, color=ARROW, lw=1.0, zorder=1)
    a, b = pts[-2], pts[-1]
    ax.annotate("", xy=b, xytext=a,
                arrowprops=dict(arrowstyle="-|>", color=ARROW, lw=1.0,
                                mutation_scale=0.7, shrinkA=0, shrinkB=0), zorder=2)


def draw_box(ax, box, title_fs, sub_fs):
    x, y, w, h = box["x"], box["y"], box["w"], box["h"]
    c = COLORS[box["kind"]]
    ax.add_patch(FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0.01,rounding_size=0.05",
        linewidth=1.1, edgecolor=c["ec"], facecolor=c["fc"], zorder=3))
    if box["sub"]:
        ax.text(x, y + 0.21 * h, box["title"], ha="center", va="center",
                fontsize=title_fs, fontweight="bold", color=c["tc"], zorder=4)
        ax.text(x, y - 0.23 * h, box["sub"], ha="center", va="center",
                fontsize=sub_fs, color="#5d6d7e", zorder=4)
    else:
        ax.text(x, y, box["title"], ha="center", va="center",
                fontsize=title_fs, fontweight="bold", color=c["tc"], zorder=4)


# edges shared by both layouts: (src, src_side, dst, dst_side)
EDGES = [
    ("input", "bottom", "check", "top"),
    ("check", "bottom", "meta", "top"),
    ("meta", "bottom", "pre", "top"),
    ("pre", "bottom", "clean", "top"),
    ("clean", "bottom", "read_based", "top"),
    ("clean", "bottom", "assembly", "top"),
    ("read_based", "bottom", "merge", "top"),
    ("assembly", "bottom", "mapping", "top"),
    ("mapping", "bottom", "binning", "top"),
    ("binning", "bottom", "mag_qc", "top"),
    ("mag_qc", "bottom", "drep", "top"),
    ("drep", "bottom", "rep", "top"),
    ("rep", "bottom", "taxonomy", "top"),
    ("rep", "bottom", "gene", "top"),
    ("rep", "bottom", "abundance", "top"),
    ("gene", "bottom", "annotation", "top"),
    ("taxonomy", "bottom", "integration", "top"),
    ("annotation", "bottom", "integration", "top"),
    ("integration", "bottom", "multiqc", "top"),
    ("integration", "bottom", "plot", "top"),
]


def draw_edges(ax, boxes, invert_y=False):
    for src, sside, dst, dside in EDGES:
        draw_arrow(ax, [anchor(boxes[src], sside, invert_y),
                        anchor(boxes[dst], dside, invert_y)])
    # abundance -> integration routed along the right margin
    ab, it = boxes["abundance"], boxes["integration"]
    draw_arrow(ax, [anchor(ab, "bottom", invert_y), (ab["x"], it["y"]),
                    anchor(it, "right", invert_y)])


# ---------------------------------------------------------------------------
# VERTICAL layout (data-unit coordinates) -> workflow.png
# ---------------------------------------------------------------------------
def render_vertical(out):
    TITLE_FS, SUB_FS = 9.0, 7.0
    H_PROC, H_DATA, X_PAD, MIN_W = 0.40, 0.30, 0.30, 1.7

    nodes = [
        ("input",      6.0, 17.3, "samplesheet.csv",       "--input",                    "data"),
        ("check",      6.0, 16.1, "CHECK_SAMPLESHEET",     None,                         "process"),
        ("meta",       6.0, 15.0, "tuple(meta, reads)",    None,                         "data"),
        ("pre",        6.0, 13.7, "PREPROCESSING",         "FASTQC → FASTP → HOST_REMOVAL", "process"),
        ("clean",      6.0, 12.6, "clean reads",           None,                         "data"),
        ("read_based", 2.5, 11.2, "READ_BASED",            "KRAKEN2 → BRACKEN ∥ HUMANN", "process"),
        ("assembly",   9.5, 11.2, "ASSEMBLY",              "MEGAHIT ∥ metaSPAdes → QUAST", "process"),
        ("merge",      2.5,  9.9, "READ_BASED_MERGE",      "sample×taxa / sample×pathway", "process"),
        ("mapping",    9.5,  9.9, "MAPPING",               "BOWTIE2 → SAMTOOLS → DEPTH",  "process"),
        ("binning",    9.5,  8.6, "BINNING",               "MetaBAT2 → SPLIT_BINS",       "process"),
        ("mag_qc",     9.5,  7.4, "MAG_QC",                "CHECKM2",                     "process"),
        ("drep",       9.5,  6.2, "DREP",                  "dRep dereplication",          "process"),
        ("rep",        9.5,  5.1, "representative MAGs",   None,                         "data"),
        ("taxonomy",   4.5,  3.8, "TAXONOMY",              "GTDB-Tk",                     "process"),
        ("gene",       7.7,  3.8, "GENE_PREDICTION",       "Prodigal",                    "process"),
        ("abundance", 10.7,  3.8, "ABUNDANCE",             "CoverM",                      "process"),
        ("annotation", 7.7,  2.7, "ANNOTATION",            "DIAMOND ∥ eggNOG ∥ RGI",      "process"),
        ("integration", 6.0, 1.55, "INTEGRATION",          "two core tables",             "output"),
        ("multiqc",     4.3, 0.55, "MULTIQC",              "summary report",              "output"),
        ("plot",        7.7, 0.55, "PLOTTING",             "7 PNG figures",               "output"),
    ]

    def text_w(text, fontsize, weight="normal"):
        t = ax.text(0, 0, text, fontsize=fontsize, fontweight=weight)
        bb = t.get_window_extent(renderer=ax.figure.canvas.get_renderer())
        t.remove()
        x0, x1 = ax.get_xlim()
        return bb.width / ax.figure.dpi * (x1 - x0) / ax.figure.get_size_inches()[0]

    fig, ax = plt.subplots(figsize=(7.8, 11.4), dpi=200)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 18)
    ax.axis("off")

    boxes = {}
    for key, x, y, title, sub, kind in nodes:
        tw = text_w(title, TITLE_FS, "bold")
        sw = text_w(sub, SUB_FS) if sub else 0
        w = max(tw, sw) + 2 * X_PAD
        w = max(w, MIN_W)
        h = (H_DATA if kind == "data" else H_PROC) * 2
        boxes[key] = dict(x=x, y=y, title=title, sub=sub, kind=kind, w=w, h=h)

    draw_edges(ax, boxes)

    ax.text(6.0, 17.7, "nextflow-metagenomics — workflow overview",
            ha="center", va="top", fontsize=11, fontweight="bold", color="#1c2833")
    ax.text(6.0, 0.04,
            "PLOTTING also reads the merged read-based matrices and upstream "
            "binning / QC / taxonomy tables.  ∥ = parallel branches.",
            ha="center", va="bottom", fontsize=6.8, color="#7f8c8d")

    for box in boxes.values():
        draw_box(ax, box, TITLE_FS, SUB_FS)

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    print(f"wrote {out}")


# ---------------------------------------------------------------------------
# COMPACT layout (inch coordinates, auto y) -> workflow_compact.png
# ---------------------------------------------------------------------------
COMPACT_LEVELS = [
    [("input",      3.50, "samplesheet.csv",       "--input",                    "data")],
    [("check",      3.50, "CHECK_SAMPLESHEET",     None,                         "process")],
    [("meta",       3.50, "tuple(meta, reads)",    None,                         "data")],
    [("pre",        3.50, "PREPROCESSING",         "FASTQC → FASTP → HOST_REMOVAL", "process")],
    [("clean",      3.50, "clean reads",           None,                         "data")],
    [("read_based", 1.35, "READ_BASED",            "KRAKEN2 → BRACKEN ∥ HUMANN", "process"),
     ("assembly",   5.65, "ASSEMBLY",              "MEGAHIT ∥ metaSPAdes → QUAST", "process")],
    [("merge",      1.35, "READ_BASED_MERGE",      "sample×taxa / sample×pathway", "process"),
     ("mapping",    5.65, "MAPPING",               "BOWTIE2 → SAMTOOLS → DEPTH",  "process")],
    [("binning",    5.65, "BINNING",               "MetaBAT2 → SPLIT_BINS",       "process")],
    [("mag_qc",     5.65, "MAG_QC",                "CHECKM2",                     "process")],
    [("drep",       5.65, "DREP",                  "dRep dereplication",          "process")],
    [("rep",        5.65, "representative MAGs",   None,                         "data")],
    [("taxonomy",   2.30, "TAXONOMY",              "GTDB-Tk",                     "process"),
     ("gene",       4.25, "GENE_PREDICTION",       "Prodigal",                    "process"),
     ("abundance",  5.90, "ABUNDANCE",             "CoverM",                      "process")],
    [("annotation", 4.25, "ANNOTATION",            "DIAMOND ∥ eggNOG ∥ RGI",      "process")],
    [("integration", 3.50, "INTEGRATION",          "two core tables",             "output")],
    [("multiqc",     2.45, "MULTIQC",              "summary report",              "output"),
     ("plot",        4.55, "PLOTTING",             "7 PNG figures",               "output")],
]


def render_compact(out):
    TITLE_FS, SUB_FS = 8.0, 6.3
    H_PROC, H_DATA, X_PAD, GAP = 0.22, 0.14, 0.24, 0.15
    FIG_W = 7.0

    def text_w_inch(ax, text, fontsize, weight="normal"):
        t = ax.text(0, 0, text, fontsize=fontsize, fontweight=weight)
        bb = t.get_window_extent(renderer=ax.figure.canvas.get_renderer())
        t.remove()
        return bb.width / ax.figure.dpi

    # temporary figure for text measurement
    tmp_fig, tmp_ax = plt.subplots(figsize=(FIG_W, 1.0))

    # measure widths + assign half-heights (in inches)
    half = {}
    widths = {}
    for level in COMPACT_LEVELS:
        for key, x, title, sub, kind in level:
            tw = text_w_inch(tmp_ax, title, TITLE_FS, "bold")
            sw = text_w_inch(tmp_ax, sub, SUB_FS) if sub else 0
            widths[key] = max(tw, sw) + 2 * X_PAD
            half[key] = H_PROC if sub else H_DATA

    plt.close(tmp_fig)

    # compute y top-down (in inches)
    boxes = {}
    top = 0.0
    for level in COMPACT_LEVELS:
        h = max(2 * half[k] for k, *_ in level)
        for key, x, title, sub, kind in level:
            boxes[key] = dict(x=x, y=top + h / 2, w=widths[key], h=2 * half[key],
                              title=title, sub=sub, kind=kind)
        top += h + GAP
    total_h = top - GAP

    fig_h = total_h + 1.0
    fig, ax = plt.subplots(figsize=(FIG_W, fig_h), dpi=200)
    ax.set_xlim(0, FIG_W)
    ax.set_ylim(total_h + 0.5, -0.5)  # inverted: y increases downward
    ax.axis("off")

    draw_edges(ax, boxes, invert_y=True)

    ax.text(FIG_W / 2, -0.25, "nextflow-metagenomics — workflow overview",
            ha="center", va="center", fontsize=10, fontweight="bold", color="#1c2833")
    ax.text(FIG_W / 2, total_h + 0.25,
            "PLOTTING also reads the merged read-based matrices and upstream "
            "binning / QC / taxonomy tables.  ∥ = parallel branches.",
            ha="center", va="center", fontsize=6.2, color="#7f8c8d")

    for box in boxes.values():
        draw_box(ax, box, TITLE_FS, SUB_FS)

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", facecolor="white")

    print(f"figsize = ({FIG_W:.2f}, {fig_h:.2f}) in")
    print("box extents (x0,x1 / y0,y1):")
    for key, b in boxes.items():
        print(f"  {key:12s} x=[{b['x']-b['w']/2:.2f},{b['x']+b['w']/2:.2f}] "
              f"y=[{b['y']-b['h']/2:.2f},{b['y']+b['h']/2:.2f}]")
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--layout", choices=["vertical", "compact"], default="vertical")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if args.layout == "vertical":
        out = args.out or REPO_ROOT / "docs/images/workflow.png"
        render_vertical(out)
    else:
        out = args.out or REPO_ROOT / "docs/images/workflow_compact.png"
        render_compact(out)


if __name__ == "__main__":
    main()
