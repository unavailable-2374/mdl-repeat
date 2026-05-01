#!/usr/bin/env python3
"""
Generate Test F: bimodal subfamily synthetic genome.

Design: ONE ancestral consensus (500 bp) gives rise to two clades:
  - Group Lo: 100 pristine copies at 2% divergence from cons_lo
  - Group Hi:  80 older copies  at 2% divergence from cons_hi
    (cons_lo and cons_hi are themselves 18% diverged from each other,
    derived from the same ancestral base)

Both groups share the same discovery seed l-mer (>80% identity to base),
so the discover stage finds them as ONE combined family.  The split stage
should detect the bimodal divergence distribution and attempt a split.

What we actually assert in the test:
  1. With -vv, "[split]" log lines are printed for the family (instrumentation works).
  2. The same number of families is produced under -threads 1 and -threads 4.
  3. No regression: total families in [1, 4] range (splits can produce more).

We do NOT assert n_splits == 1 because MDL may correctly reject the split
if the pre-check (using old num_edits) shows insufficient gain — which is
by design (MDL is the correct arbiter).

Output: tests/data/testF.fa
Usage:  python3 tests/gen_testF.py -o tests/data
"""

import random
import sys
import os
import argparse


def random_seq(length, rng):
    return ''.join(rng.choice(list('ACGT')) for _ in range(length))


def mutate_seq(seq, divergence, rng):
    bases = list(seq)
    n_mut = int(len(bases) * divergence)
    positions = rng.sample(range(len(bases)), min(n_mut, len(bases)))
    for pos in positions:
        orig = bases[pos]
        bases[pos] = rng.choice([b for b in 'ACGT' if b != orig])
    return ''.join(bases)


def write_fasta(filename, name, sequence, line_width=80):
    with open(filename, 'w') as f:
        f.write(f'>{name}\n')
        for i in range(0, len(sequence), line_width):
            f.write(sequence[i:i+line_width] + '\n')


def write_bed(filename, annotations):
    with open(filename, 'w') as f:
        for chrom, start, end, name, strand in annotations:
            f.write(f'{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n')


def gen_testF(output_dir):
    rng = random.Random(99)
    genome_len = 1_000_000
    chrom = 'testF'

    # Two completely independent consensi that share 80% identity
    # (enough to be merged by the merge stage into one family)
    # Then discovered together, they form one bimodal family.
    cons_base = random_seq(500, rng)

    # Lo subfamily: 2% divergence from base (very clean)
    cons_lo = mutate_seq(cons_base, 0.02, rng)
    # Hi subfamily: 18% divergence from base (clearly older)
    cons_hi = mutate_seq(cons_base, 0.18, rng)

    genome = list(random_seq(genome_len, rng))
    annotations = []
    used_ranges = []

    def find_free(length, tries=3000):
        for _ in range(tries):
            pos = rng.randint(0, genome_len - length - 1)
            ok = all(not (pos < e and pos + length > s) for (s, e) in used_ranges)
            if ok:
                used_ranges.append((pos, pos + length))
                return pos
        return None

    # Group Lo: 25 copies at 2% divergence from cons_lo (= ~4% from cons_base)
    lo_count = 0
    for _ in range(25):
        pos = find_free(500)
        if pos is None:
            continue
        copy = mutate_seq(cons_lo, 0.02, rng)
        genome[pos:pos + 500] = list(copy)
        annotations.append((chrom, pos, pos + 500, 'Fam_Lo', '+'))
        lo_count += 1

    # Group Hi: 20 copies at 2% divergence from cons_hi (= ~20% from cons_base)
    hi_count = 0
    for _ in range(20):
        pos = find_free(500)
        if pos is None:
            continue
        copy = mutate_seq(cons_hi, 0.02, rng)
        genome[pos:pos + 500] = list(copy)
        annotations.append((chrom, pos, pos + 500, 'Fam_Hi', '+'))
        hi_count += 1

    genome_str = ''.join(genome)
    write_fasta(os.path.join(output_dir, 'testF.fa'), chrom, genome_str)
    write_bed(os.path.join(output_dir, 'testF_truth.bed'),
              sorted(annotations, key=lambda x: x[1]))
    print(f"Test F: {len(genome_str)} bp genome, "
          f"Lo={lo_count} copies (4% from base), Hi={hi_count} copies (20% from base)")
    return lo_count, hi_count


def main():
    parser = argparse.ArgumentParser(description='Generate testF bimodal-subfamily genome')
    parser.add_argument('-o', '--output', default='tests/data')
    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)
    gen_testF(args.output)


if __name__ == '__main__':
    main()
