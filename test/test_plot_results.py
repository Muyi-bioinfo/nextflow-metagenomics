#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_plot_results.py — bin/plot_results.py 单测 (Phase 21)

纯标准库 unittest (需 matplotlib + numpy), 直接运行:
    python3 test/test_plot_results.py

覆盖场景 (对应 STATUS 测试矩阵): PNG 非空 + 尺寸 (IHDR 800×600) / 缺列报错 /
空表单点不崩 / PCoA 缺失跳过 / 漏斗缺失层级 / function top-N。合成表仅用于
验证"代码能画" —— 真实 database-dependent 数据本机无库, 不伪造 (规则 3)。
"""
import os
import struct
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import plot_results as pr  # noqa: E402


def _write(path, content):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(content)
    return path


def _png_dims(path):
    """读 PNG IHDR 的宽高 (像素)。"""
    with open(path, "rb") as fh:
        assert fh.read(8) == b"\x89PNG\r\n\x1a\n"
        assert fh.read(4) == b"\x00\x00\x00\r"  # IHDR length
        assert fh.read(4) == b"IHDR"
        w = struct.unpack(">I", fh.read(4))[0]
        h = struct.unpack(">I", fh.read(4))[0]
    return w, h


def _run(argv):
    """跑主入口 (parse_args + func), 返回 (rc, 是否产出 PNG)。"""
    rc = pr.main(argv)
    return rc, os.path.exists(argv[argv.index("--output") + 1])


MAG_ABUNDANCE = "MAG_ID\tS01\tS02\nS02.megahit.metabat2.001\t0\t0.57479572\n"
QC_HEADER = "meta_id\tmag_id\tName\tCompleteness\tContamination\n"
TAX_HEADER = ("sample\tmag_id\trep_mag_id\tdomain\tphylum\tclass\torder\t"
              "family\tgenus\tspecies\n")
BIN_HEADER = ("mag_id\tassembly_unit\tassembler\tassembly_mode\tbinner\t"
              "source_bin\tn_contigs\ttotal_bp\tlargest_contig_bp\t"
              "mean_contig_bp\tgc_percent\tn_bases\n")


class TestAbundanceHeatmap(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_produces_png_with_dimensions(self):
        p = _write(os.path.join(self.tmp, "ab.tsv"), MAG_ABUNDANCE)
        out = os.path.join(self.tmp, "f.png")
        rc, ok = _run(["mag-abundance-heatmap", "--matrix", p,
                       "--output", out])
        self.assertEqual(rc, 0)
        self.assertTrue(ok)
        self.assertGreater(os.path.getsize(out), 0)
        self.assertEqual(_png_dims(out), (800, 600))

    def test_empty_table_skips_no_png(self):
        p = _write(os.path.join(self.tmp, "ab.tsv"), "")  # 0 字节
        out = os.path.join(self.tmp, "f.png")
        rc, ok = _run(["mag-abundance-heatmap", "--matrix", p,
                       "--output", out])
        self.assertEqual(rc, 0)
        self.assertFalse(ok)  # 跳过, 不写 PNG

    def test_header_only_skips(self):
        p = _write(os.path.join(self.tmp, "ab.tsv"), "MAG_ID\tS01\tS02\n")
        out = os.path.join(self.tmp, "f.png")
        rc, ok = _run(["mag-abundance-heatmap", "--matrix", p,
                       "--output", out])
        self.assertEqual(rc, 0)
        self.assertFalse(ok)


class TestQcScatter(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_missing_column_raises(self):
        p = _write(os.path.join(self.tmp, "qc.tsv"),
                   "meta_id\tmag_id\tName\nS02\tm1\tx\n")
        out = os.path.join(self.tmp, "f.png")
        with self.assertRaises(SystemExit):
            pr.main(["qc-scatter", "--qc", p, "--output", out])

    def test_single_point_does_not_crash(self):
        p = _write(os.path.join(self.tmp, "qc.tsv"),
                   QC_HEADER + "S02\tm1\tx\t75.0\t5.0\n")
        out = os.path.join(self.tmp, "f.png")
        rc, ok = _run(["qc-scatter", "--qc", p, "--output", out])
        self.assertEqual(rc, 0)
        self.assertTrue(ok)

    def test_empty_table_skips(self):
        p = _write(os.path.join(self.tmp, "qc.tsv"), QC_HEADER)  # 仅表头
        out = os.path.join(self.tmp, "f.png")
        rc, ok = _run(["qc-scatter", "--qc", p, "--output", out])
        self.assertEqual(rc, 0)
        self.assertFalse(ok)


class TestTaxonomyComposition(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_missing_phylum_column_raises(self):
        p = _write(os.path.join(self.tmp, "tax.tsv"),
                   "sample\tmag_id\tdomain\nS02\tm1\tBacteria\n")
        out = os.path.join(self.tmp, "f.png")
        with self.assertRaises(SystemExit):
            pr.main(["mag-taxonomy-composition", "--taxonomy", p,
                     "--output", out])

    def test_composition_draws(self):
        p = _write(os.path.join(self.tmp, "tax.tsv"),
                   TAX_HEADER +
                   "S02\tm1\t\tBacteria\tFirmicutes\tBacilli\t\t\t\t\n"
                   "S01\tm2\tm1\tBacteria\tFirmicutes\tBacilli\t\t\t\t\n")
        out = os.path.join(self.tmp, "f.png")
        rc, ok = _run(["mag-taxonomy-composition", "--taxonomy", p,
                       "--output", out])
        self.assertEqual(rc, 0)
        self.assertTrue(ok)


class TestTaxonomicComposition(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_top_n_via_helper(self):
        # 直接测 top_labels: top-N + Other 合并
        labels, values = pr.top_labels(["a", "b", "c", "d"], [10, 5, 3, 2], 2)
        self.assertEqual(labels, ["a", "b", "Other"])
        self.assertEqual(values, [10, 5, 5])

    def test_picks_s_level(self):
        # merged_S + merged_G 同时给, 应选 S 层级
        ps = _write(os.path.join(self.tmp, "merged_S.tsv"),
                    "name\tS01\tS02\nA\t0.5\t0.2\n")
        _write(os.path.join(self.tmp, "merged_G.tsv"),
               "name\tS01\tS02\nB\t0.3\t0.1\n")
        out = os.path.join(self.tmp, "f.png")
        rc, ok = _run(["taxonomic-composition", "--matrices",
                       ps, os.path.join(self.tmp, "merged_G.tsv"),
                       "--output", out])
        self.assertEqual(rc, 0)
        self.assertTrue(ok)


class TestBetaPcoa(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_missing_distance_skips(self):
        out = os.path.join(self.tmp, "f.png")
        rc, ok = _run(["beta-pcoa", "--distance",
                       os.path.join(self.tmp, "nope.tsv"), "--output", out])
        self.assertEqual(rc, 0)
        self.assertFalse(ok)

    def test_valid_pcoa_draws(self):
        p = _write(os.path.join(self.tmp, "beta.tsv"),
                   "\tS01\tS02\tS03\n"
                   "S01\t0\t0.5\t1.0\n"
                   "S02\t0.5\t0\t0.6\n"
                   "S03\t1.0\t0.6\t0\n")
        out = os.path.join(self.tmp, "f.png")
        rc, ok = _run(["beta-pcoa", "--distance", p, "--output", out])
        self.assertEqual(rc, 0)
        self.assertTrue(ok)

    def test_two_samples_skip(self):
        # 2 样本 Bray-Curtis 只有 1 个正特征值 → 跳过
        p = _write(os.path.join(self.tmp, "beta.tsv"),
                   "\tS01\tS02\nS01\t0\t0.5\nS02\t0.5\t0\n")
        out = os.path.join(self.tmp, "f.png")
        rc, ok = _run(["beta-pcoa", "--distance", p, "--output", out])
        self.assertEqual(rc, 0)
        self.assertFalse(ok)


class TestPathwayHeatmap(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_top_n_draws(self):
        body = "".join(f"PWY-{i}: pathway {i}\t{10 - i}\t{i}\n"
                       for i in range(5))
        p = _write(os.path.join(self.tmp, "pw.tsv"),
                   "pathway\tS01\tS02\n" + body)
        out = os.path.join(self.tmp, "f.png")
        rc, ok = _run(["pathway-heatmap", "--matrix", p, "--output", out,
                       "--top", "3"])
        self.assertEqual(rc, 0)
        self.assertTrue(ok)

    def test_empty_skips(self):
        p = _write(os.path.join(self.tmp, "pw.tsv"), "pathway\tS01\tS02\n")
        out = os.path.join(self.tmp, "f.png")
        rc, ok = _run(["pathway-heatmap", "--matrix", p, "--output", out])
        self.assertEqual(rc, 0)
        self.assertFalse(ok)


class TestWorkflowSummary(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def _paths(self, with_membership=True, with_taxonomy=True):
        binp = _write(os.path.join(self.tmp, "bin.tsv"),
                      BIN_HEADER + "m1\tS01\tmegahit\tsingle\tmetabat2\tx.fa\t"
                      "10\t50000\t8000\t5000.0\t45.2\t0\n")
        qcp = _write(os.path.join(self.tmp, "qc.tsv"),
                     QC_HEADER + "S02\tm1\tx\t75.0\t5.0\n")
        mem = _write(os.path.join(self.tmp, "mem.tsv"),
                     "m1\tS01\tm1\nm2\tS02\tm1\n") if with_membership else ""
        tax = _write(os.path.join(self.tmp, "tax.tsv"),
                     TAX_HEADER +
                     "S02\tm1\t\tBacteria\tFirmicutes\tBacilli\t\t\t\t\n"
                     ) if with_taxonomy else ""
        return binp, qcp, mem, tax

    def test_all_levels_draws(self):
        b, q, m, t = self._paths()
        out = os.path.join(self.tmp, "f.png")
        rc, ok = _run(["workflow-summary", "--bin-summary", b, "--qc", q,
                       "--membership", m, "--taxonomy", t, "--output", out])
        self.assertEqual(rc, 0)
        self.assertTrue(ok)

    def test_missing_middle_level_skips(self):
        # 无 membership (0 字节) → "After dRep" 层级跳过, 其余画
        b, q, _, t = self._paths(with_membership=False)
        mem = _write(os.path.join(self.tmp, "mem_empty.tsv"), "")
        out = os.path.join(self.tmp, "f.png")
        rc, ok = _run(["workflow-summary", "--bin-summary", b, "--qc", q,
                       "--membership", mem, "--taxonomy", t, "--output", out])
        self.assertEqual(rc, 0)
        self.assertTrue(ok)

    def test_all_missing_skips(self):
        out = os.path.join(self.tmp, "f.png")
        rc, ok = _run(["workflow-summary",
                       "--bin-summary", os.path.join(self.tmp, "nb.tsv"),
                       "--qc", os.path.join(self.tmp, "nq.tsv"),
                       "--membership", os.path.join(self.tmp, "nm.tsv"),
                       "--taxonomy", os.path.join(self.tmp, "nt.tsv"),
                       "--output", out])
        self.assertEqual(rc, 0)
        self.assertFalse(ok)


class TestLoadMatrix(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_returns_none_for_missing(self):
        self.assertIsNone(pr.load_matrix(os.path.join(self.tmp, "nope.tsv")))

    def test_returns_none_for_header_only(self):
        p = _write(os.path.join(self.tmp, "m.tsv"), "a\tS01\tS02\n")
        self.assertIsNone(pr.load_matrix(p))

    def test_parses_values(self):
        p = _write(os.path.join(self.tmp, "m.tsv"),
                   "a\tS01\tS02\nx\t0.5\t0\n")
        labels, cols, data = pr.load_matrix(p)
        self.assertEqual(labels, ["x"])
        self.assertEqual(cols, ["S01", "S02"])
        self.assertEqual(data.tolist(), [[0.5, 0.0]])


if __name__ == "__main__":
    unittest.main(verbosity=2)
