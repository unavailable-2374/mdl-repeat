#!/usr/bin/env python3
"""Bounded rescue/refine parameter tuning harness for mdl-repeat fixtures.

This script runs mdl-repeat on committed test fixtures and scores predicted
instance BED files against the fixture truth BED files. Metrics printed by this
script are fixture/test metrics only; they are not biological benchmark claims.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

TOOLS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TOOLS_DIR.parent
sys.path.insert(0, str(TOOLS_DIR))

from bed_pr import BedInterval, evaluate, load_bed  # noqa: E402


@dataclass(frozen=True)
class Fixture:
    name: str
    fasta: Path
    truth_bed: Path


@dataclass(frozen=True)
class TuningCase:
    name: str
    args: tuple[str, ...]


@dataclass(frozen=True)
class BaseMetrics:
    truth_bp: int
    predicted_bp: int
    tp_bp: int

    @property
    def precision(self) -> float:
        return self.tp_bp / self.predicted_bp if self.predicted_bp else 0.0

    @property
    def recall(self) -> float:
        return self.tp_bp / self.truth_bp if self.truth_bp else 0.0

    @property
    def f1(self) -> float:
        denom = self.precision + self.recall
        return 2.0 * self.precision * self.recall / denom if denom else 0.0


FIXTURES: dict[str, Fixture] = {
    "testA": Fixture(
        "testA", REPO_ROOT / "tests/data/testA.fa", REPO_ROOT / "tests/data/testA_truth.bed"
    ),
    "testB": Fixture(
        "testB", REPO_ROOT / "tests/data/testB.fa", REPO_ROOT / "tests/data/testB_truth.bed"
    ),
    "testC": Fixture(
        "testC", REPO_ROOT / "tests/data/testC.fa", REPO_ROOT / "tests/data/testC_truth.bed"
    ),
    "testD": Fixture(
        "testD", REPO_ROOT / "tests/data/testD.fa", REPO_ROOT / "tests/data/testD_truth.bed"
    ),
    "testF": Fixture(
        "testF", REPO_ROOT / "tests/data/testF.fa", REPO_ROOT / "tests/data/testF_truth.bed"
    ),
    "testG": Fixture(
        "testG", REPO_ROOT / "tests/data/testG.fa", REPO_ROOT / "tests/data/testG_truth.bed"
    ),
    "multichr": Fixture(
        "multichr",
        REPO_ROOT / "tests/data/multichr/multichr.fa",
        REPO_ROOT / "tests/data/multichr/multichr_truth.bed",
    ),
}

PRESETS: dict[str, tuple[TuningCase, ...]] = {
    "smoke": (
        TuningCase("baseline", ()),
        TuningCase(
            "rescue_bounded",
            (
                "-recall-rescue",
                "-rescue-l-delta",
                "1",
                "-rescue-maxrepeats",
                "5",
                "-rescue-min-gap",
                "100",
            ),
        ),
        TuningCase(
            "refine_relaxed",
            (
                "-max-divergence",
                "0.35",
                "-refine-gap",
                "-3",
                "-refine-maxoffset",
                "16",
            ),
        ),
        TuningCase(
            "rescue_refine_relaxed",
            (
                "-recall-rescue",
                "-rescue-l-delta",
                "1",
                "-rescue-maxrepeats",
                "5",
                "-rescue-min-gap",
                "100",
                "-max-divergence",
                "0.35",
                "-refine-gap",
                "-3",
                "-refine-maxoffset",
                "16",
            ),
        ),
    ),
    "small": (
        TuningCase("baseline", ()),
        TuningCase(
            "rescue_l_delta_1",
            (
                "-recall-rescue",
                "-rescue-l-delta",
                "1",
                "-rescue-maxrepeats",
                "5",
                "-rescue-min-gap",
                "100",
            ),
        ),
        TuningCase(
            "rescue_l_delta_2",
            (
                "-recall-rescue",
                "-rescue-l-delta",
                "2",
                "-rescue-maxrepeats",
                "5",
                "-rescue-min-gap",
                "100",
            ),
        ),
        TuningCase("max_divergence_025", ("-max-divergence", "0.25")),
        TuningCase("max_divergence_035", ("-max-divergence", "0.35")),
        TuningCase("refine_gap_minus3", ("-refine-gap", "-3")),
        TuningCase("refine_maxoffset_16", ("-refine-maxoffset", "16")),
        TuningCase("max_dp_cells_5000000", ("-max-dp-cells", "5000000")),
        TuningCase("coalesce_factor_0", ("-coalesce-factor", "0")),
        TuningCase("coalesce_factor_10", ("-coalesce-factor", "10")),
    ),
}

TSV_FIELDS = (
    "metric_scope",
    "fixture",
    "case",
    "status",
    "exit_code",
    "elapsed_sec",
    "timed_out",
    "families",
    "instances",
    "truth_instances",
    "interval_tp",
    "interval_fp",
    "interval_fn",
    "interval_precision",
    "interval_recall",
    "interval_f1",
    "truth_bp",
    "predicted_bp",
    "tp_bp",
    "bp_precision",
    "bp_recall",
    "bp_f1",
    "command",
    "stdout_log",
    "stderr_log",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run bounded mdl-repeat rescue/refine parameter sweeps on committed fixtures. "
            "All reported metrics are fixture/test metrics."
        )
    )
    parser.add_argument(
        "--binary",
        type=Path,
        default=REPO_ROOT / "bin/mdl-repeat",
        help="Path to mdl-repeat binary (default: bin/mdl-repeat)",
    )
    parser.add_argument(
        "--preset",
        choices=sorted(PRESETS),
        default="smoke",
        help="Tuning matrix preset. Default is deliberately small.",
    )
    parser.add_argument(
        "--fixtures",
        default="testC,multichr",
        help=(
            "Comma-separated fixtures to run, or 'all'. "
            f"Available: {','.join(sorted(FIXTURES))}. Default: testC,multichr"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "tests/results/tune_rescue_refine",
        help="Directory for per-run outputs and summary TSV.",
    )
    parser.add_argument(
        "--timeout-sec",
        type=float,
        default=60.0,
        help="Per mdl-repeat run timeout in seconds (default: 60).",
    )
    parser.add_argument(
        "--min-overlap",
        type=float,
        default=0.5,
        help="Interval overlap threshold for BED instance metrics (default: 0.5).",
    )
    parser.add_argument(
        "--overlap-mode",
        choices=("pred", "truth", "min"),
        default="min",
        help="Denominator for interval overlap fraction (default: min).",
    )
    parser.add_argument(
        "--ignore-strand",
        action="store_true",
        help="Ignore strand when scoring interval and bp-level fixture metrics.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Pass -threads N to mdl-repeat for each run (default: 1).",
    )
    parser.add_argument(
        "--limit-runs",
        type=int,
        default=None,
        help="Optional hard cap on total runs, useful for local smoke checks.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned commands without running mdl-repeat.",
    )
    return parser.parse_args()


def selected_fixtures(names_arg: str) -> list[Fixture]:
    if names_arg == "all":
        return [FIXTURES[name] for name in sorted(FIXTURES)]
    names = [name.strip() for name in names_arg.split(",") if name.strip()]
    unknown = [name for name in names if name not in FIXTURES]
    if unknown:
        raise ValueError(f"unknown fixture(s): {','.join(unknown)}")
    return [FIXTURES[name] for name in names]


def validate_inputs(binary: Path, fixtures: list[Fixture]) -> None:
    if not binary.exists() or not binary.is_file():
        raise FileNotFoundError(
            f"mdl-repeat binary not found: {binary}. Build it first with `make`."
        )
    missing: list[Path] = []
    for fixture in fixtures:
        if not fixture.fasta.exists():
            missing.append(fixture.fasta)
        if not fixture.truth_bed.exists():
            missing.append(fixture.truth_bed)
    if missing:
        joined = "\n  ".join(str(path) for path in missing)
        raise FileNotFoundError(f"fixture input file(s) missing:\n  {joined}")


def count_fasta_records(path: Path) -> int:
    if not path.exists():
        return 0
    count = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def interval_key(interval: BedInterval, ignore_strand: bool) -> tuple[str, str]:
    return (interval.chrom, "*") if ignore_strand else (interval.chrom, interval.strand)


def merged_segments(intervals: list[BedInterval], ignore_strand: bool) -> dict[tuple[str, str], list[tuple[int, int]]]:
    buckets: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for interval in intervals:
        buckets.setdefault(interval_key(interval, ignore_strand), []).append(
            (interval.start, interval.end)
        )

    merged: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for key, segments in buckets.items():
        segments.sort()
        out: list[tuple[int, int]] = []
        for start, end in segments:
            if not out or start > out[-1][1]:
                out.append((start, end))
            else:
                old_start, old_end = out[-1]
                out[-1] = (old_start, max(old_end, end))
        merged[key] = out
    return merged


def segment_length(segments_by_key: dict[tuple[str, str], list[tuple[int, int]]]) -> int:
    return sum(end - start for segments in segments_by_key.values() for start, end in segments)


def intersect_segment_lists(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> int:
    i = 0
    j = 0
    total = 0
    while i < len(a) and j < len(b):
        start = max(a[i][0], b[j][0])
        end = min(a[i][1], b[j][1])
        if end > start:
            total += end - start
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return total


def base_metrics(
    predicted: list[BedInterval], truth: list[BedInterval], ignore_strand: bool
) -> BaseMetrics:
    pred_merged = merged_segments(predicted, ignore_strand)
    truth_merged = merged_segments(truth, ignore_strand)
    truth_bp = segment_length(truth_merged)
    predicted_bp = segment_length(pred_merged)
    tp_bp = 0
    for key, pred_segments in pred_merged.items():
        truth_segments = truth_merged.get(key)
        if truth_segments:
            tp_bp += intersect_segment_lists(pred_segments, truth_segments)
    return BaseMetrics(truth_bp=truth_bp, predicted_bp=predicted_bp, tp_bp=tp_bp)


def fmt_float(value: float) -> str:
    return f"{value:.4f}"


def command_for_run(
    binary: Path,
    fixture: Fixture,
    case: TuningCase,
    out_fa: Path,
    out_bed: Path,
    stats_tsv: Path,
    stderr_log: Path,
    threads: int,
) -> list[str]:
    cmd = [
        str(binary),
        "-sequence",
        str(fixture.fasta),
        "-output",
        str(out_fa),
        "-instances",
        str(out_bed),
        "-stats",
        str(stats_tsv),
        "-threads",
        str(threads),
    ]
    cmd.extend(case.args)
    if "-recall-rescue" in case.args:
        audit_tsv = stderr_log.with_suffix(".rescue_audit.tsv")
        cmd.extend(["-rescue-audit", str(audit_tsv)])
    return cmd


def empty_row(
    fixture: Fixture,
    case: TuningCase,
    status: str,
    exit_code: int | str,
    elapsed_sec: float,
    timed_out: bool,
    cmd: list[str],
    stdout_log: Path,
    stderr_log: Path,
) -> dict[str, str]:
    row = {field: "" for field in TSV_FIELDS}
    row.update(
        {
            "metric_scope": "fixture/test",
            "fixture": fixture.name,
            "case": case.name,
            "status": status,
            "exit_code": str(exit_code),
            "elapsed_sec": fmt_float(elapsed_sec),
            "timed_out": "true" if timed_out else "false",
            "command": " ".join(cmd),
            "stdout_log": str(stdout_log),
            "stderr_log": str(stderr_log),
        }
    )
    return row


def run_case(
    binary: Path,
    fixture: Fixture,
    case: TuningCase,
    output_dir: Path,
    timeout_sec: float,
    min_overlap: float,
    overlap_mode: str,
    ignore_strand: bool,
    threads: int,
    dry_run: bool,
) -> dict[str, str]:
    run_dir = output_dir / fixture.name / case.name
    run_dir.mkdir(parents=True, exist_ok=True)
    out_fa = run_dir / "families.fa"
    out_bed = run_dir / "instances.bed"
    stats_tsv = run_dir / "stats.tsv"
    stdout_log = run_dir / "stdout.log"
    stderr_log = run_dir / "stderr.log"
    cmd = command_for_run(
        binary, fixture, case, out_fa, out_bed, stats_tsv, stderr_log, threads
    )

    if dry_run:
        row = empty_row(
            fixture, case, "dry-run", "", 0.0, False, cmd, stdout_log, stderr_log
        )
        print("DRY-RUN\t" + row["command"])
        return row

    start = time.monotonic()
    timed_out = False
    try:
        completed = subprocess.run(
            cmd,
            cwd=REPO_ROOT,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout_sec,
        )
        exit_code: int | str = completed.returncode
        stdout_log.write_text(completed.stdout)
        stderr_log.write_text(completed.stderr)
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        exit_code = "timeout"
        stdout_log.write_text(exc.stdout or "")
        stderr_log.write_text(exc.stderr or "")
    elapsed_sec = time.monotonic() - start

    if timed_out or exit_code != 0:
        return empty_row(
            fixture,
            case,
            "timeout" if timed_out else "failed",
            exit_code,
            elapsed_sec,
            timed_out,
            cmd,
            stdout_log,
            stderr_log,
        )

    predicted = load_bed(out_bed) if out_bed.exists() else []
    truth = load_bed(fixture.truth_bed)
    interval, _ = evaluate(
        predicted,
        truth,
        min_overlap=min_overlap,
        overlap_mode=overlap_mode,
        ignore_strand=ignore_strand,
        by_class=False,
    )
    bases = base_metrics(predicted, truth, ignore_strand)

    row = empty_row(
        fixture,
        case,
        "ok",
        exit_code,
        elapsed_sec,
        timed_out,
        cmd,
        stdout_log,
        stderr_log,
    )
    row.update(
        {
            "families": str(count_fasta_records(out_fa)),
            "instances": str(len(predicted)),
            "truth_instances": str(len(truth)),
            "interval_tp": str(interval.tp),
            "interval_fp": str(interval.fp),
            "interval_fn": str(interval.fn),
            "interval_precision": fmt_float(interval.precision),
            "interval_recall": fmt_float(interval.recall),
            "interval_f1": fmt_float(interval.f1),
            "truth_bp": str(bases.truth_bp),
            "predicted_bp": str(bases.predicted_bp),
            "tp_bp": str(bases.tp_bp),
            "bp_precision": fmt_float(bases.precision),
            "bp_recall": fmt_float(bases.recall),
            "bp_f1": fmt_float(bases.f1),
        }
    )
    return row


def print_row(row: dict[str, str]) -> None:
    fields = [
        row["metric_scope"],
        row["fixture"],
        row["case"],
        row["status"],
        f"families={row['families'] or 'NA'}",
        f"instances={row['instances'] or 'NA'}",
        f"interval_f1={row['interval_f1'] or 'NA'}",
        f"bp_f1={row['bp_f1'] or 'NA'}",
        f"elapsed={row['elapsed_sec']}s",
    ]
    print("\t".join(fields), flush=True)


def main() -> int:
    args = parse_args()
    try:
        fixtures = selected_fixtures(args.fixtures)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    cases = list(PRESETS[args.preset])
    if args.limit_runs is not None:
        cases_by_fixture = [(fixture, case) for fixture in fixtures for case in cases]
        cases_by_fixture = cases_by_fixture[: args.limit_runs]
    else:
        cases_by_fixture = [(fixture, case) for fixture in fixtures for case in cases]

    try:
        validate_inputs(args.binary, fixtures)
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_tsv = args.output_dir / "summary.tsv"

    print(
        "NOTE: all values below are fixture/test metrics computed from committed "
        "tests/data truth BED files; they are not biological benchmark results.",
        flush=True,
    )
    print(
        f"Preset: {args.preset}; fixtures: {','.join(f.name for f in fixtures)}",
        flush=True,
    )
    print(f"Output directory: {args.output_dir}", flush=True)

    rows: list[dict[str, str]] = []
    for fixture, case in cases_by_fixture:
        row = run_case(
            args.binary,
            fixture,
            case,
            args.output_dir,
            args.timeout_sec,
            args.min_overlap,
            args.overlap_mode,
            args.ignore_strand,
            args.threads,
            args.dry_run,
        )
        rows.append(row)
        print_row(row)

    with summary_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote fixture/test metric table: {summary_tsv}", flush=True)

    failed = [row for row in rows if row["status"] not in {"ok", "dry-run"}]
    if failed:
        print(
            f"ERROR: {len(failed)} run(s) failed or timed out; inspect stderr_log paths.",
            flush=True,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
