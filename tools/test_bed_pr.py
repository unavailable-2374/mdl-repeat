#!/usr/bin/env python3
"""
Unit tests for tools/bed_pr.py.

Run:  python3 tools/test_bed_pr.py
Exit: 0 = all pass, 1 = at least one assertion failed.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

# Make `bed_pr` importable when run from repo root or tools/ dir
sys.path.insert(0, str(Path(__file__).resolve().parent))

from bed_pr import (BedInterval, BedIndex, overlap_fraction, evaluate,
                    load_bed, main)

n_pass = 0
n_fail = 0


def check(cond: bool, msg: str) -> None:
    global n_pass, n_fail
    if cond:
        print(f"  PASS: {msg}")
        n_pass += 1
    else:
        print(f"  FAIL: {msg}")
        n_fail += 1


# ---------------------------------------------------------------------------

def test_overlap_fraction():
    print("\n[1] overlap_fraction")
    a = BedInterval("chr1", 100, 200, "p1", "+")
    b = BedInterval("chr1", 150, 250, "t1", "+")
    # Overlap is 50; pred len 100; truth len 100; min len 100
    check(abs(overlap_fraction(a, b, "pred") - 0.5) < 1e-9, "pred-mode = 0.5")
    check(abs(overlap_fraction(a, b, "truth") - 0.5) < 1e-9, "truth-mode = 0.5")
    check(abs(overlap_fraction(a, b, "min") - 0.5) < 1e-9, "min-mode = 0.5")
    # Disjoint
    c = BedInterval("chr1", 300, 400, "t2", "+")
    check(overlap_fraction(a, c, "min") == 0.0, "disjoint = 0")


def test_index():
    print("\n[2] BedIndex strand-aware lookup")
    truth = [
        BedInterval("chr1", 100, 200, "L1",   "+"),
        BedInterval("chr1", 150, 250, "L1",   "-"),
        BedInterval("chr1", 1000, 2000, "L2", "+"),
    ]
    idx = BedIndex(truth, ignore_strand=False)
    q_plus  = BedInterval("chr1", 120, 180, "p", "+")
    hits_p  = idx.overlaps(q_plus)
    check(len(hits_p) == 1 and hits_p[0].name == "L1" and hits_p[0].strand == "+",
          "strand-aware: + query hits + truth only")

    idx_any = BedIndex(truth, ignore_strand=True)
    hits_any = idx_any.overlaps(q_plus)
    check(len(hits_any) == 2, "ignore-strand: query hits both + and - truth")


def test_perfect_match():
    print("\n[3] perfect match → P=R=F1=1")
    pred = [
        BedInterval("chr1", 100, 200, "F1", "+"),
        BedInterval("chr1", 1000, 1100, "F1", "+"),
    ]
    truth = [
        BedInterval("chr1", 100, 200, "FAM", "+"),
        BedInterval("chr1", 1000, 1100, "FAM", "+"),
    ]
    overall, _ = evaluate(pred, truth, min_overlap=0.5,
                          overlap_mode="min", ignore_strand=False,
                          by_class=False)
    check(overall.tp == 2 and overall.fp == 0 and overall.fn == 0,
          f"TP=2 FP=0 FN=0 (got TP={overall.tp} FP={overall.fp} FN={overall.fn})")
    check(overall.precision == 1.0 and overall.recall == 1.0,
          "precision = recall = 1.0")


def test_partial_overlap_threshold():
    print("\n[4] overlap below threshold counts as FP+FN")
    # 30% overlap — fails default min-overlap=0.5
    pred = [BedInterval("chr1", 100, 200, "F1", "+")]
    truth = [BedInterval("chr1", 170, 250, "FAM", "+")]
    overall, _ = evaluate(pred, truth, min_overlap=0.5, overlap_mode="min",
                          ignore_strand=False, by_class=False)
    check(overall.tp == 0 and overall.fp == 1 and overall.fn == 1,
          f"low-overlap → FP+FN (got TP={overall.tp} FP={overall.fp} "
          f"FN={overall.fn})")
    # Lower threshold and it becomes a TP
    overall2, _ = evaluate(pred, truth, min_overlap=0.2, overlap_mode="min",
                           ignore_strand=False, by_class=False)
    check(overall2.tp == 1 and overall2.fp == 0 and overall2.fn == 0,
          "with min-overlap=0.2 it becomes TP")


def test_by_class():
    print("\n[5] per-class breakdown")
    pred = [
        BedInterval("chr1", 100, 200, "p", "+"),  # match L1
        BedInterval("chr1", 500, 600, "p", "+"),  # spurious
    ]
    truth = [
        BedInterval("chr1", 100, 200, "L1", "+"),
        BedInterval("chr1", 800, 900, "L2", "+"),  # missed
    ]
    overall, by_cls = evaluate(pred, truth, min_overlap=0.5,
                               overlap_mode="min", ignore_strand=False,
                               by_class=True)
    check(overall.tp == 1 and overall.fp == 1 and overall.fn == 1,
          f"overall TP=1 FP=1 FN=1 (got {overall.tp}/{overall.fp}/{overall.fn})")
    assert by_cls is not None
    check(by_cls["L1"].tp == 1 and by_cls["L1"].fn == 0,
          "L1: TP=1 FN=0 (matched)")
    check(by_cls["L2"].tp == 0 and by_cls["L2"].fn == 1,
          "L2: TP=0 FN=1 (missed)")


def test_load_bed_skip_invalid():
    print("\n[6] load_bed skips comments / malformed rows")
    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as fh:
        fh.write("# header comment\n")
        fh.write("track name=foo\n")
        fh.write("chr1\t100\t200\tF1\t0\t+\n")
        fh.write("chr1\t300\tINVALID\tF2\t0\t+\n")
        fh.write("\n")
        fh.write("chr1\t500\t400\tF3\t0\t+\n")  # end < start
        fh.write("chr1\t700\t800\tF4\t0\t-\n")
        path = Path(fh.name)
    try:
        items = load_bed(path)
        check(len(items) == 2, f"loaded 2 valid rows (got {len(items)})")
        check(items[0].chrom == "chr1" and items[0].start == 100,
              "first row parsed correctly")
        check(items[1].strand == "-", "strand parsed correctly")
    finally:
        path.unlink()


def test_main_against_synth_truth():
    print("\n[7] CLI smoke against tests/data/testA_truth.bed if present")
    truth_path = Path(__file__).resolve().parents[1] / "tests/data/testA_truth.bed"
    if not truth_path.exists():
        check(True, "(skipped — testA_truth.bed not generated)")
        return
    # Compare truth against itself: must give P=R=F1=1.0
    rc = main(["--predicted", str(truth_path),
               "--truth",     str(truth_path),
               "--min-overlap", "0.5",
               "--output", "/dev/null"])
    check(rc == 0, "main() returns 0 on valid inputs")


# ---------------------------------------------------------------------------

def main_test() -> int:
    print("=== bed_pr.py unit tests ===")
    test_overlap_fraction()
    test_index()
    test_perfect_match()
    test_partial_overlap_threshold()
    test_by_class()
    test_load_bed_skip_invalid()
    test_main_against_synth_truth()
    print(f"\n{n_pass} passed, {n_fail} failed")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main_test())
