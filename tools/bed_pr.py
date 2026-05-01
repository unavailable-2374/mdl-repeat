#!/usr/bin/env python3
"""
bed_pr.py — Compute precision / recall / F1 of predicted repeat instances
against a ground-truth BED file.

Designed for benchmarking mdl-repeat output against:
  - synthetic ground truth (tests/data/test*_truth.bed)
  - Dfam annotations on hg38 (e.g. chr19) — see README

A predicted instance is counted as TP for a truth class iff there exists
a truth interval of that class on the same chromosome+strand whose
overlap fraction (per the user-chosen mode: relative to predicted, truth,
or smaller) exceeds --min-overlap.  Interval indexing uses sorted-sweep
in O((n + m) log n) — no external libraries.

Output is a tab-separated summary, one row per truth class plus an
overall row, suitable for piping into a spreadsheet or further analysis.

Usage:
    bed_pr.py --predicted predicted.bed --truth truth.bed \\
              [--min-overlap 0.5] [--overlap-mode pred|truth|min] \\
              [--ignore-strand] [--by-class] \\
              [--output report.tsv]

Exit codes:
    0  = success
    2  = invalid input
"""

from __future__ import annotations

import argparse
import bisect
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class BedInterval:
    chrom: str
    start: int        # 0-based inclusive
    end: int          # 0-based exclusive
    name: str
    strand: str       # '+', '-', or '.'

    @property
    def length(self) -> int:
        return self.end - self.start


def load_bed(path: Path) -> list[BedInterval]:
    """Read a BED6-or-better file.  Truncates extra columns silently.
    Skips browser/track lines and empty / comment lines."""
    out: list[BedInterval] = []
    with path.open() as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track") \
                    or line.startswith("browser"):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                sys.stderr.write(f"{path}:{lineno}: <3 columns, skipping\n")
                continue
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError:
                sys.stderr.write(f"{path}:{lineno}: non-integer coords, "
                                 "skipping\n")
                continue
            if end <= start:
                sys.stderr.write(f"{path}:{lineno}: end <= start, skipping\n")
                continue
            name   = cols[3] if len(cols) >= 4 else "."
            strand = cols[5] if len(cols) >= 6 else "."
            if strand not in ("+", "-", "."):
                strand = "."
            out.append(BedInterval(cols[0], start, end, name, strand))
    return out


# ---------------------------------------------------------------------------
# Index for fast overlap lookup: for each (chrom, strand) bucket, store
# intervals sorted by start.  An overlap query at [s, e) finds the
# first interval whose start >= e (upper bound), then scans backward as
# long as the interval's end > s.
# ---------------------------------------------------------------------------

class BedIndex:
    def __init__(self, intervals: list[BedInterval], ignore_strand: bool):
        self._buckets: dict[tuple[str, str], list[BedInterval]] = \
            defaultdict(list)
        self._starts: dict[tuple[str, str], list[int]] = defaultdict(list)
        self._ignore_strand = ignore_strand
        for it in intervals:
            key = self._key(it.chrom, it.strand)
            self._buckets[key].append(it)
        for key in self._buckets:
            self._buckets[key].sort(key=lambda x: (x.start, x.end))
            self._starts[key] = [x.start for x in self._buckets[key]]
            # Maintain a parallel running-max of end so we can stop scans
            # early.  (Simple O(n) precompute.)
        # We also need max_end_seen across leftward scans, but for clarity
        # we just bound the scan at min_start - max_truth_len.
        self._max_len: dict[tuple[str, str], int] = {}
        for key, lst in self._buckets.items():
            self._max_len[key] = max((it.length for it in lst), default=0)

    def _key(self, chrom: str, strand: str) -> tuple[str, str]:
        return (chrom, "*") if self._ignore_strand else (chrom, strand)

    def overlaps(self, query: BedInterval) -> list[BedInterval]:
        """Return all intervals overlapping the query (open-end semantics)."""
        keys = []
        if self._ignore_strand:
            keys.append((query.chrom, "*"))
        else:
            keys.append((query.chrom, query.strand))
            if query.strand == ".":
                # query strandless: try both
                keys.extend([(query.chrom, "+"), (query.chrom, "-")])
        out: list[BedInterval] = []
        for key in keys:
            bucket = self._buckets.get(key)
            if not bucket:
                continue
            starts = self._starts[key]
            # Right boundary: first start >= query.end
            hi = bisect.bisect_left(starts, query.end)
            # Left boundary: any interval with start >= query.start - max_len
            lo_bound = query.start - self._max_len[key]
            lo = bisect.bisect_left(starts, lo_bound)
            for i in range(lo, hi):
                it = bucket[i]
                if it.end > query.start:
                    out.append(it)
        return out


# ---------------------------------------------------------------------------
# P/R computation
# ---------------------------------------------------------------------------

def overlap_fraction(a: BedInterval, b: BedInterval, mode: str) -> float:
    """Fraction of overlap relative to the chosen denominator."""
    s = max(a.start, b.start)
    e = min(a.end,   b.end)
    ov = max(0, e - s)
    if mode == "pred":
        denom = a.length
    elif mode == "truth":
        denom = b.length
    elif mode == "min":
        denom = min(a.length, b.length)
    else:
        raise ValueError(f"unknown overlap mode: {mode}")
    return ov / denom if denom > 0 else 0.0


@dataclass
class ClassMetrics:
    tp: int = 0
    fp: int = 0
    fn: int = 0

    @property
    def precision(self) -> float:
        denom = self.tp + self.fp
        return self.tp / denom if denom > 0 else 0.0

    @property
    def recall(self) -> float:
        denom = self.tp + self.fn
        return self.tp / denom if denom > 0 else 0.0

    @property
    def f1(self) -> float:
        p, r = self.precision, self.recall
        return (2 * p * r / (p + r)) if (p + r) > 0 else 0.0


def evaluate(predicted: list[BedInterval],
             truth: list[BedInterval],
             min_overlap: float,
             overlap_mode: str,
             ignore_strand: bool,
             by_class: bool):
    truth_idx = BedIndex(truth, ignore_strand)
    pred_idx  = BedIndex(predicted, ignore_strand)

    # Per-class buckets — class derived from truth name.  The "predicted"
    # class is whatever the predicted BED's column 4 says, but for class-
    # specific metrics we can only score predictions that match a truth
    # class via overlap (or we treat all predictions as a single bucket).
    metrics: dict[str, ClassMetrics] = defaultdict(ClassMetrics)
    overall = ClassMetrics()

    # TP / FP from the predicted side
    for p in predicted:
        hits = truth_idx.overlaps(p)
        matched = False
        matched_classes: set[str] = set()
        for t in hits:
            if overlap_fraction(p, t, overlap_mode) >= min_overlap:
                matched = True
                matched_classes.add(t.name)
        if matched:
            overall.tp += 1
            for cls in matched_classes:
                metrics[cls].tp += 1
        else:
            overall.fp += 1
            metrics["__unmatched_pred__"].fp += 1

    # FN from the truth side
    for t in truth:
        hits = pred_idx.overlaps(t)
        matched = any(
            overlap_fraction(p, t, overlap_mode) >= min_overlap
            for p in hits
        )
        if not matched:
            overall.fn += 1
            metrics[t.name].fn += 1

    return overall, metrics if by_class else None


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def fmt_row(name: str, m: ClassMetrics) -> str:
    return (f"{name}\t{m.tp}\t{m.fp}\t{m.fn}\t"
            f"{m.precision:.4f}\t{m.recall:.4f}\t{m.f1:.4f}")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--predicted", required=True, type=Path)
    p.add_argument("--truth",     required=True, type=Path)
    p.add_argument("--min-overlap", type=float, default=0.5,
                   help="Required reciprocal-overlap fraction (default 0.5)")
    p.add_argument("--overlap-mode", choices=("pred", "truth", "min"),
                   default="min",
                   help="Denominator for overlap fraction (default: min)")
    p.add_argument("--ignore-strand", action="store_true",
                   help="Ignore strand when matching")
    p.add_argument("--by-class", action="store_true",
                   help="Emit per-truth-class breakdown in addition to overall")
    p.add_argument("--output", type=Path, default=None,
                   help="Write TSV report to this file (default: stdout)")
    args = p.parse_args(argv)

    if not 0.0 < args.min_overlap <= 1.0:
        sys.stderr.write("--min-overlap must be in (0, 1]\n")
        return 2

    if not args.predicted.exists():
        sys.stderr.write(f"predicted not found: {args.predicted}\n")
        return 2
    if not args.truth.exists():
        sys.stderr.write(f"truth not found: {args.truth}\n")
        return 2

    predicted = load_bed(args.predicted)
    truth     = load_bed(args.truth)

    overall, by_cls = evaluate(predicted, truth,
                               args.min_overlap, args.overlap_mode,
                               args.ignore_strand, args.by_class)

    out_fp = args.output.open("w") if args.output else sys.stdout
    try:
        out_fp.write("class\tTP\tFP\tFN\tprecision\trecall\tF1\n")
        out_fp.write(fmt_row("__overall__", overall) + "\n")
        if by_cls is not None:
            for cls in sorted(by_cls.keys()):
                out_fp.write(fmt_row(cls, by_cls[cls]) + "\n")
    finally:
        if args.output:
            out_fp.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
