#!/usr/bin/env python3
"""
Generate a 2-chromosome synthetic FASTA for ENG-N8/9/11 cross-chromosome
correctness tests.

Layout:
  chr1 (100 kb): random background + 6 copies of REPEAT_A (1 kb each),
                 6 copies of REPEAT_B (800 bp), spaced uniformly.  The
                 last REPEAT_A copy is placed ~500 bp before chr1 end —
                 close enough to chr2's first repeat to trip the ENG-N8
                 tandem-coalesce / ENG-N9 fragment-assembly bug if the
                 cross-chromosome guard is missing.

  chr2 (100 kb): same scheme but with REPEAT_C (1 kb) and REPEAT_D (800 bp).
                 First REPEAT_C copy is placed ~500 bp into chr2 — close
                 enough to chr1's tail repeat for cross-boundary mishaps.

Random seed is fixed so the file is reproducible.
"""
import random
import sys
import os

random.seed(20260428)

BASES = "ACGT"

def randseq(n):
    return "".join(random.choice(BASES) for _ in range(n))

def mutate(seq, divergence):
    s = list(seq)
    n = len(s)
    n_mut = int(n * divergence)
    for _ in range(n_mut):
        i = random.randrange(n)
        s[i] = random.choice([b for b in BASES if b != s[i]])
    return "".join(s)

def revcomp(s):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp[b] for b in reversed(s))

# --- Repeat consensus sequences (distinct enough to be separate families)
REPEAT_A = randseq(1000)   # chr1 family A
REPEAT_B = randseq(800)    # chr1 family B
REPEAT_C = randseq(1000)   # chr2 family C
REPEAT_D = randseq(800)    # chr2 family D

# --- Build chr1 -------------------------------------------------------
CHR1_LEN  = 100_000
CHR2_LEN  = 100_000

def build_chrom(length, repeat_specs, seed_offset):
    """
    repeat_specs: list of (consensus, n_copies, divergence, position_hints)
      position_hints: list of integer start positions; len must equal n_copies
    Returns (sequence_string, list of (name, start, end, strand))
    """
    rng = random.Random(seed_offset)
    # backbone: random
    backbone = list("".join(rng.choice(BASES) for _ in range(length)))
    intervals = []  # list of (name, start, end, strand)
    for name, cons, n_copies, div, positions, strands in repeat_specs:
        for i, pos in enumerate(positions):
            strand = strands[i % len(strands)]
            inst = mutate(cons, div)
            if strand == "-":
                inst = revcomp(inst)
            end = pos + len(inst)
            if end > length:
                continue
            backbone[pos:end] = list(inst)
            intervals.append((name, pos, end, strand))
    return "".join(backbone), intervals


# chr1 layout: spread 6 copies of A and 6 copies of B uniformly
# Last A copy ~500 bp before chr1 end (position 98_500-99_500)
chr1_A_pos = [5_000, 20_000, 38_000, 55_000, 72_000, 98_500]
chr1_B_pos = [12_000, 28_000, 45_000, 62_000, 80_000, 92_000]

chr1_seq, chr1_iv = build_chrom(
    CHR1_LEN,
    [
        ("RepA", REPEAT_A, 6, 0.05, chr1_A_pos, ["+", "+", "-", "+", "-", "+"]),
        ("RepB", REPEAT_B, 6, 0.05, chr1_B_pos, ["+", "-", "+", "+", "-", "-"]),
    ],
    seed_offset=1,
)

# chr2 layout: 6 copies of C, 6 copies of D
# First C copy at position 500 — ~ near chr2 start, to test cross-boundary
chr2_C_pos = [500, 18_000, 35_000, 52_000, 70_000, 89_000]
chr2_D_pos = [10_000, 25_000, 42_000, 60_000, 78_000, 95_000]

chr2_seq, chr2_iv = build_chrom(
    CHR2_LEN,
    [
        ("RepC", REPEAT_C, 6, 0.05, chr2_C_pos, ["+", "-", "+", "+", "+", "-"]),
        ("RepD", REPEAT_D, 6, 0.05, chr2_D_pos, ["-", "+", "+", "-", "+", "+"]),
    ],
    seed_offset=2,
)

# --- Write FASTA -----------------------------------------------------
out_dir = os.path.dirname(os.path.abspath(__file__))
fa_path = os.path.join(out_dir, "multichr.fa")
truth_path = os.path.join(out_dir, "multichr_truth.bed")

with open(fa_path, "w") as fh:
    fh.write(">chr1\n")
    for i in range(0, len(chr1_seq), 60):
        fh.write(chr1_seq[i:i+60] + "\n")
    fh.write(">chr2\n")
    for i in range(0, len(chr2_seq), 60):
        fh.write(chr2_seq[i:i+60] + "\n")

with open(truth_path, "w") as fh:
    for name, s, e, strand in chr1_iv:
        fh.write(f"chr1\t{s}\t{e}\t{name}\t.\t{strand}\n")
    for name, s, e, strand in chr2_iv:
        fh.write(f"chr2\t{s}\t{e}\t{name}\t.\t{strand}\n")

print(f"Wrote {fa_path}: chr1={len(chr1_seq)} bp, chr2={len(chr2_seq)} bp")
print(f"Wrote {truth_path}: {len(chr1_iv) + len(chr2_iv)} truth intervals")
