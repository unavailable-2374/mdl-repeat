#!/usr/bin/env python3
"""
Generate synthetic genomes with ground-truth repeat annotations.

Creates FASTA + ground-truth BED files for testing mdl-repeat.
"""

import random
import sys
import os
import argparse

def random_seq(length):
    """Generate a random DNA sequence."""
    return ''.join(random.choice('ACGT') for _ in range(length))

def mutate_seq(seq, divergence):
    """Mutate a sequence with given divergence fraction."""
    bases = list(seq)
    n_mut = int(len(bases) * divergence)
    positions = random.sample(range(len(bases)), min(n_mut, len(bases)))
    for pos in positions:
        orig = bases[pos]
        bases[pos] = random.choice([b for b in 'ACGT' if b != orig])
    return ''.join(bases)

def write_fasta(filename, name, sequence, line_width=80):
    """Write a single-sequence FASTA file."""
    with open(filename, 'w') as f:
        f.write(f'>{name}\n')
        for i in range(0, len(sequence), line_width):
            f.write(sequence[i:i+line_width] + '\n')

def write_bed(filename, annotations):
    """Write BED file with annotations. Each entry: (chrom, start, end, name, strand)."""
    with open(filename, 'w') as f:
        for chrom, start, end, name, strand in annotations:
            f.write(f'{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n')


def test_a(output_dir):
    """Test A: 1MB genome, 3 families.
    - Family 1: 50 copies of 300bp element at 0% divergence
    - Family 2: 20 copies of 500bp element at 10% divergence
    - Family 3: 100 copies of 150bp element at 5% divergence
    """
    random.seed(42)
    genome_len = 1000000
    chrom = 'testA'

    # Generate consensus sequences
    fam1_cons = random_seq(300)
    fam2_cons = random_seq(500)
    fam3_cons = random_seq(150)

    # Generate background
    genome = list(random_seq(genome_len))
    annotations = []

    # Insert copies at random positions (non-overlapping)
    used_ranges = []

    def find_free_pos(length, tries=1000):
        for _ in range(tries):
            pos = random.randint(0, genome_len - length - 1)
            overlap = False
            for (s, e) in used_ranges:
                if pos < e and pos + length > s:
                    overlap = True
                    break
            if not overlap:
                used_ranges.append((pos, pos + length))
                return pos
        return None

    # Family 1: 50 copies, 0% divergence
    for i in range(50):
        pos = find_free_pos(300)
        if pos is None: continue
        copy = fam1_cons  # 0% divergence
        genome[pos:pos+300] = list(copy)
        annotations.append((chrom, pos, pos+300, 'Family1', '+'))

    # Family 2: 20 copies, 10% divergence
    for i in range(20):
        pos = find_free_pos(500)
        if pos is None: continue
        copy = mutate_seq(fam2_cons, 0.10)
        genome[pos:pos+500] = list(copy)
        annotations.append((chrom, pos, pos+500, 'Family2', '+'))

    # Family 3: 100 copies, 5% divergence
    for i in range(100):
        pos = find_free_pos(150)
        if pos is None: continue
        copy = mutate_seq(fam3_cons, 0.05)
        genome[pos:pos+150] = list(copy)
        annotations.append((chrom, pos, pos+150, 'Family3', '+'))

    genome_str = ''.join(genome)
    write_fasta(os.path.join(output_dir, 'testA.fa'), chrom, genome_str)
    write_bed(os.path.join(output_dir, 'testA_truth.bed'),
              sorted(annotations, key=lambda x: x[1]))
    print(f"Test A: {len(genome_str)} bp genome, {len(annotations)} repeat instances")


def test_b(output_dir):
    """Test B: Tandem array.
    - 50bp unit repeated 15 times with 3% divergence
    - Embedded in 100KB background
    """
    random.seed(43)
    genome_len = 100000
    chrom = 'testB'

    unit = random_seq(50)
    tandem_copies = 15
    tandem_div = 0.03

    genome = list(random_seq(genome_len))
    annotations = []

    # Insert tandem array at position 50000
    pos = 50000
    for i in range(tandem_copies):
        copy = mutate_seq(unit, tandem_div)
        start = pos + i * 50
        genome[start:start+50] = list(copy)
        annotations.append((chrom, start, start+50, 'TandemUnit', '+'))

    genome_str = ''.join(genome)
    write_fasta(os.path.join(output_dir, 'testB.fa'), chrom, genome_str)
    write_bed(os.path.join(output_dir, 'testB_truth.bed'),
              sorted(annotations, key=lambda x: x[1]))
    print(f"Test B: {len(genome_str)} bp genome, tandem array of {tandem_copies}x{50}bp")


def test_c(output_dir):
    """Test C: Near detection limit.
    - 3 copies of 200bp at 0% divergence (should be detected)
    - 2 copies of 200bp at 20% divergence (should be borderline/rejected)
    """
    random.seed(44)
    genome_len = 200000
    chrom = 'testC'

    cons_good = random_seq(200)
    cons_bad = random_seq(200)

    genome = list(random_seq(genome_len))
    annotations = []

    # 3 copies at 0% divergence
    positions = [10000, 50000, 90000]
    for pos in positions:
        genome[pos:pos+200] = list(cons_good)
        annotations.append((chrom, pos, pos+200, 'GoodFamily', '+'))

    # 2 copies at 20% divergence
    positions2 = [130000, 170000]
    for pos in positions2:
        copy = mutate_seq(cons_bad, 0.20)
        genome[pos:pos+200] = list(copy)
        annotations.append((chrom, pos, pos+200, 'BorderlineFamily', '+'))

    genome_str = ''.join(genome)
    write_fasta(os.path.join(output_dir, 'testC.fa'), chrom, genome_str)
    write_bed(os.path.join(output_dir, 'testC_truth.bed'),
              sorted(annotations, key=lambda x: x[1]))
    print(f"Test C: {len(genome_str)} bp genome, detection limit test")


def test_d(output_dir):
    """Test D: Nested TE (SINE inside LINE).
    - LINE: 1000bp element, 10 copies at 8% divergence
    - SINE: 200bp element, inserted into each LINE at position 400-600
    - Total SINE copies: 10 (inside LINEs) + 20 solo copies
    """
    random.seed(45)
    genome_len = 500000
    chrom = 'testD'

    line_cons = random_seq(1000)
    sine_cons = random_seq(200)

    # Create nested LINE (SINE inside)
    nested_line = line_cons[:400] + sine_cons + line_cons[600:]

    genome = list(random_seq(genome_len))
    annotations = []
    used_ranges = []

    def find_free(length, tries=1000):
        for _ in range(tries):
            pos = random.randint(0, genome_len - length - 1)
            ok = True
            for (s, e) in used_ranges:
                if pos < e and pos + length > s:
                    ok = False
                    break
            if ok:
                used_ranges.append((pos, pos + length))
                return pos
        return None

    # 10 copies of nested LINE (contains SINE)
    for i in range(10):
        pos = find_free(1000)
        if pos is None: continue
        copy = mutate_seq(nested_line, 0.08)
        genome[pos:pos+1000] = list(copy)
        annotations.append((chrom, pos, pos+1000, 'LINE', '+'))
        annotations.append((chrom, pos+400, pos+600, 'SINE_nested', '+'))

    # 20 solo copies of SINE
    for i in range(20):
        pos = find_free(200)
        if pos is None: continue
        copy = mutate_seq(sine_cons, 0.05)
        genome[pos:pos+200] = list(copy)
        annotations.append((chrom, pos, pos+200, 'SINE_solo', '+'))

    genome_str = ''.join(genome)
    write_fasta(os.path.join(output_dir, 'testD.fa'), chrom, genome_str)
    write_bed(os.path.join(output_dir, 'testD_truth.bed'),
              sorted(annotations, key=lambda x: x[1]))
    print(f"Test D: {len(genome_str)} bp genome, nested TE test")


def main():
    parser = argparse.ArgumentParser(description='Generate synthetic test genomes')
    parser.add_argument('-o', '--output', default='tests/data',
                        help='Output directory (default: tests/data)')
    parser.add_argument('-t', '--test', choices=['A', 'B', 'C', 'D', 'all'],
                        default='all', help='Which test to generate')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    if args.test in ('A', 'all'): test_a(args.output)
    if args.test in ('B', 'all'): test_b(args.output)
    if args.test in ('C', 'all'): test_c(args.output)
    if args.test in ('D', 'all'): test_d(args.output)

    print(f"\nAll test data written to {args.output}/")


if __name__ == '__main__':
    main()
