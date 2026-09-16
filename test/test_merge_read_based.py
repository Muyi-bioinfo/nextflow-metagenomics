#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_merge_read_based.py — bin/merge_read_based.py 单测 (Phase 20)

纯标准库 unittest, 直接运行:  python3 test/test_merge_read_based.py
覆盖场景 空表 / 缺列 / 缺样本补 0 / 层级选择 /
单样本距离矩阵跳过 / 数值与逐样本一致。合成表仅用于验证合并逻辑 —— 真实
Bracken/HUMAnN 数据本地无库, 不虚构 。
"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import merge_read_based as mrb  # noqa: E402


def _write(path, content):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(content)
    return path


def _read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return [l.rstrip("\n").split("\t") for l in fh]


BRACKEN_HEADER = "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\t" \
                 "added_reads\tnew_est_reads\tfraction_total_reads\n"


class TestParseBracken(unittest.TestCase):
    def test_parse_normal(self):
        p = _write(os.path.join(self._tmp, "s.bracken.S.tsv"),
                   BRACKEN_HEADER +
                   "Streptococcus\t1311\tS\t100\t0\t100\t0.50000\n"
                   "Escherichia\t561\tS\t50\t0\t50\t0.25000\n")
        rows = mrb.parse_bracken(p)
        self.assertEqual(rows, [("Streptococcus", "0.50000"),
                                ("Escherichia", "0.25000")])

    def test_parse_empty_table(self):
        # 0 字节表 (stub touch / 0-read 样本) → 空列表, 不报错
        p = _write(os.path.join(self._tmp, "empty.tsv"), "")
        self.assertEqual(mrb.parse_bracken(p), [])

    def test_parse_header_only(self):
        p = _write(os.path.join(self._tmp, "hdr.tsv"), BRACKEN_HEADER)
        self.assertEqual(mrb.parse_bracken(p), [])

    def test_parse_missing_column(self):
        p = _write(os.path.join(self._tmp, "bad.tsv"),
                   "name\ttaxonomy_id\nA\t1\n")
        with self.assertRaises(SystemExit):
            mrb.parse_bracken(p)

    def setUp(self):
        self._tmp = tempfile.mkdtemp()


class TestParsePathabundance(unittest.TestCase):
    def test_parse_skips_comment_header(self):
        p = _write(os.path.join(self._tmp, "s_pathabundance.tsv"),
                   "# Pathway\tAbundance\n"
                   "PWY-3781: aerobic respiration I\t123.45\n"
                   "UNMAPPED\t0.0\n")
        self.assertEqual(mrb.parse_pathabundance(p),
                         [("PWY-3781: aerobic respiration I", "123.45"),
                          ("UNMAPPED", "0.0")])

    def test_parse_empty(self):
        p = _write(os.path.join(self._tmp, "empty.tsv"), "")
        self.assertEqual(mrb.parse_pathabundance(p), [])

    def test_parse_missing_abundance_col(self):
        p = _write(os.path.join(self._tmp, "bad.tsv"),
                   "# Pathway\tAbundance\nPWY-1\n")
        with self.assertRaises(SystemExit):
            mrb.parse_pathabundance(p)

    def setUp(self):
        self._tmp = tempfile.mkdtemp()


class TestMergeBracken(unittest.TestCase):
    """端到端: 走 cmd_bracken 主流程 (合成多样本/多层级表)。"""

    def _run(self, levels, samples, values):
        """构建合成 bracken 文件并跑 cmd_bracken。
        values: {(level, sample): [(name, frac), ...]}"""
        manifest_lines = []
        for level in levels:
            for sample in samples:
                fname = f"{sample}.bracken.{level}.tsv"
                rows = values.get((level, sample), [])
                body = BRACKEN_HEADER + "".join(
                    f"{name}\t0\t{level}\t0\t0\t0\t{frac}\n" for name, frac in rows)
                _write(os.path.join(self._tmp, fname), body)
                manifest_lines.append(f"{level}\t{sample}\t{os.path.join(self._tmp, fname)}\n")
        manifest = _write(os.path.join(self._tmp, "manifest.tsv"),
                          "".join(manifest_lines))
        outdir = os.path.join(self._tmp, "out")
        os.makedirs(outdir, exist_ok=True)
        rc = mrb.cmd_bracken(type("A", (), {
            "manifest": manifest, "output_dir": outdir})())
        return outdir, rc

    def test_values_match_per_sample_and_missing_filled_zero(self):
        # 数值与逐样本一致 + 缺样本补 0
        outdir, rc = self._run(
            ["S"], ["S01", "S02"],
            {("S", "S01"): [("A", "0.50000")],
             ("S", "S02"): [("A", "0.75000"), ("B", "0.25000")]})
        self.assertEqual(rc, 0)
        rows = _read_tsv(os.path.join(outdir, "merged_S.tsv"))
        self.assertEqual(rows[0], ["name", "S01", "S02"])
        data = {r[0]: r[1:] for r in rows[1:]}
        self.assertEqual(data["A"], ["0.50000", "0.75000"])  # 逐样本原样
        self.assertEqual(data["B"], ["0", "0.25000"])        # S01 缺 B → 补 0

    def test_level_selection_two_levels(self):
        outdir, rc = self._run(
            ["S", "G"], ["S01", "S02"],
            {("S", "S01"): [("SpeciesA", "0.5")],
             ("G", "S01"): [("GenusA", "0.6")]})
        self.assertEqual(rc, 0)
        self.assertTrue(os.path.exists(os.path.join(outdir, "merged_S.tsv")))
        self.assertTrue(os.path.exists(os.path.join(outdir, "merged_G.tsv")))
        # beta diversity 选 S 层级 (S 前缀优先)
        self.assertTrue(os.path.exists(os.path.join(outdir, "beta_diversity.tsv")))

    def test_beta_diversity_values(self):
        # S01 = S02 完全相同 → 距离 0; 与 S03 完全不同 → 距离 1
        outdir, rc = self._run(
            ["S"], ["S01", "S02", "S03"],
            {("S", "S01"): [("A", "0.5"), ("B", "0.5")],
             ("S", "S02"): [("A", "0.5"), ("B", "0.5")],
             ("S", "S03"): [("C", "1.0")]})
        self.assertEqual(rc, 0)
        rows = _read_tsv(os.path.join(outdir, "beta_diversity.tsv"))
        self.assertEqual(rows[0], ["", "S01", "S02", "S03"])
        idx = {rows[i][0]: i for i in range(1, len(rows))}
        d = [[float(c) for c in r[1:]] for r in rows[1:]]
        self.assertAlmostEqual(d[idx["S01"] - 1][idx["S02"] - 1], 0.0, places=6)
        self.assertAlmostEqual(d[idx["S01"] - 1][idx["S03"] - 1], 1.0, places=6)

    def test_single_sample_skips_beta(self):
        outdir, rc = self._run(
            ["S"], ["S01"],
            {("S", "S01"): [("A", "0.5"), ("B", "0.5")]})
        self.assertEqual(rc, 0)
        self.assertTrue(os.path.exists(os.path.join(outdir, "merged_S.tsv")))
        self.assertFalse(os.path.exists(os.path.join(outdir, "beta_diversity.tsv")))

    def test_empty_tables_header_only_matrix(self):
        # 全部 0 字节表 → 仅表头宽表, 不崩溃 (stub 语义)
        outdir, rc = self._run(["S"], ["S01", "S02"], {})
        self.assertEqual(rc, 0)
        rows = _read_tsv(os.path.join(outdir, "merged_S.tsv"))
        self.assertEqual(rows, [["name", "S01", "S02"]])

    def setUp(self):
        self._tmp = tempfile.mkdtemp()


class TestMergePathabundance(unittest.TestCase):
    def test_merge(self):
        tmp = tempfile.mkdtemp()
        p1 = _write(os.path.join(tmp, "S01_pathabundance.tsv"),
                    "# Pathway\tAbundance\nPWY-A: alpha\t10.0\nPWY-B: beta\t5.0\n")
        p2 = _write(os.path.join(tmp, "S02_pathabundance.tsv"),
                    "# Pathway\tAbundance\nPWY-A: alpha\t20.0\n")
        manifest = _write(os.path.join(tmp, "m.tsv"),
                          f"S01\t{p1}\nS02\t{p2}\n")
        out = os.path.join(tmp, "merged_pathabundance.tsv")
        rc = mrb.cmd_pathabundance(type("A", (), {
            "manifest": manifest, "output": out})())
        self.assertEqual(rc, 0)
        rows = _read_tsv(out)
        self.assertEqual(rows[0], ["pathway", "S01", "S02"])
        data = {r[0]: r[1:] for r in rows[1:]}
        self.assertEqual(data["PWY-A: alpha"], ["10.0", "20.0"])
        self.assertEqual(data["PWY-B: beta"], ["5.0", "0"])  # S02 缺 → 补 0


if __name__ == "__main__":
    unittest.main(verbosity=2)
