#!/usr/bin/env python3
"""
Generate Test G: 250 bp short element, 5 copies with 5% indel divergence.
Tests BLAST-based instance recruitment for short families (Phase 5 Track 1 / F').

Truth: 5 copies of a 250 bp element in a 300 kb genome.
  - 2 copies on + strand
  - 2 copies on - strand (inserted as reverse-complement)
  - 1 additional + strand copy
  Divergence: ~5% substitutions + indels.

Usage: python3 gen_testG.py -o <output_dir>
"""

import random
import os
import argparse


def random_seq(n):
    return ''.join(random.choice('ACGT') for _ in range(n))


def rc(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))


def mutate_indel(seq, div=0.05):
    """Apply substitutions + indels at given total divergence rate."""
    bases = list(seq)
    n_events = max(1, int(len(bases) * div))
    positions = sorted(
        random.sample(range(len(bases)), min(n_events, len(bases))),
        reverse=True
    )
    for pos in positions:
        op = random.choice(['sub', 'ins', 'del'])
        if op == 'sub':
            orig = bases[pos]
            bases[pos] = random.choice([b for b in 'ACGT' if b != orig])
        elif op == 'ins':
            bases.insert(pos, random.choice('ACGT'))
        else:  # del
            if len(bases) > 1:
                bases.pop(pos)
    return ''.join(bases)


def write_fasta(path, name, seq, width=80):
    with open(path, 'w') as f:
        f.write(f'>{name}\n')
        for i in range(0, len(seq), width):
            f.write(seq[i:i+width] + '\n')


def write_bed(path, records):
    with open(path, 'w') as f:
        for chrom, start, end, name, score, strand in records:
            f.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n')


def gen_testG(output_dir):
    random.seed(99)
    genome_len = 300000
    chrom = 'testG'
    cons_len = 250

    cons = random_seq(cons_len)
    genome = list(random_seq(genome_len))
    annotations = []
    used_ranges = []

    def find_free(length, tries=1000):
        for _ in range(tries):
            pos = random.randint(0, genome_len - length - 1)
            ok = all(not (pos < e and pos + length > s) for (s, e) in used_ranges)
            if ok:
                used_ranges.append((pos, pos + length))
                return pos
        return None

    # 5 copies: 3 forward (+), 2 reverse-complement (-)
    strand_choices = ['+', '+', '-', '-', '+']
    for strand in strand_choices:
        copy = mutate_indel(cons, 0.05)
        copy_len = len(copy)
        pos = find_free(copy_len)
        if pos is None:
            continue
        seq_to_insert = rc(copy) if strand == '-' else copy
        genome[pos:pos+copy_len] = list(seq_to_insert)
        annotations.append((chrom, pos, pos+copy_len, 'ShortElem', 0, strand))

    genome_str = ''.join(genome)

    write_fasta(os.path.join(output_dir, 'testG.fa'), chrom, genome_str)
    write_bed(os.path.join(output_dir, 'testG_truth.bed'),
              sorted(annotations, key=lambda x: x[1]))
    # Also write the original consensus for manual verification
    write_fasta(os.path.join(output_dir, 'testG_cons.fa'),
                f'consensus len={cons_len}', cons)

    n = len(annotations)
    print(f"Test G: {len(genome_str)} bp genome, {n} copies of {cons_len}bp element (5% indel div)")
    for (c, s, e, nm, sc, st) in sorted(annotations, key=lambda x: x[1]):
        print(f"  {st} strand: {s}-{e}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate Test G synthetic data')
    parser.add_argument('-o', '--output', default='.', help='Output directory')
    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)
    gen_testG(args.output)
