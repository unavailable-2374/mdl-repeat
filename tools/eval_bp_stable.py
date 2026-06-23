#!/usr/bin/env python3
"""Stable bp-level evaluator for TE library quality.

Bypasses RepeatMasker's global cross-match dedup, which we discovered creates
artifactual mask drops when the library is modified. Each library family is
independently rmblastn'd against the genome; hits are merged into a single BED;
recall/precision/F1 are computed against a truth BED.

Usage:
  python3 tools/eval_bp_stable.py <library.fa> <genome.fa> <truth.bed> <label>
"""
import os, sys, subprocess

if len(sys.argv) < 5:
    print(__doc__, file=sys.stderr)
    sys.exit(1)

LIB_FA = sys.argv[1]
GENOME = sys.argv[2]
TRUTH_BED = sys.argv[3]
LABEL = sys.argv[4]

RMBLASTN = os.environ.get('RMBLASTN', 'rmblastn')
MAKEBLASTDB = os.environ.get('MAKEBLASTDB', 'makeblastdb')
BEDTOOLS = os.environ.get('BEDTOOLS', 'bedtools')
WORKDIR = os.environ.get('EVAL_WORKDIR', f'/tmp/eval_stable_{LABEL}')
os.makedirs(WORKDIR, exist_ok=True)

db_path = os.path.join(WORKDIR, 'genome_db')
if not os.path.exists(db_path + '.nhr'):
    subprocess.run([MAKEBLASTDB, '-in', GENOME, '-dbtype', 'nucl',
                    '-out', db_path, '-parse_seqids'],
                   check=True, capture_output=True)

blast_out = os.path.join(WORKDIR, 'lib_vs_genome.tsv')
if not os.path.exists(blast_out):
    print(f"[{LABEL}] rmblastn ...", file=sys.stderr)
    subprocess.run([
        RMBLASTN, '-query', LIB_FA, '-db', db_path,
        '-outfmt', '6 qseqid sseqid pident length qlen slen sstart send qstart qend evalue',
        '-evalue', '1e-10', '-num_threads', '4',
        '-perc_identity', '80',
        '-gapopen', '5', '-gapextend', '2',
        '-dust', 'no', '-soft_masking', 'false',
        '-out', blast_out,
    ], check=True, capture_output=True)

mask_bed = os.path.join(WORKDIR, 'mask.bed')
with open(blast_out) as fin, open(mask_bed, 'w') as fout:
    for ln in fin:
        p = ln.rstrip().split('\t')
        if len(p) < 11: continue
        sid = p[1]
        sstart, send = int(p[6]), int(p[7])
        s, e = (sstart, send) if sstart < send else (send, sstart)
        fout.write(f"{sid}\t{s-1}\t{e}\n")

sorted_bed = os.path.join(WORKDIR, 'mask_sorted.bed')
merged_bed = os.path.join(WORKDIR, 'mask_merged.bed')
subprocess.run(f"sort -k1,1 -k2,2n {mask_bed} > {sorted_bed}", shell=True, check=True)
subprocess.run([BEDTOOLS, 'merge', '-i', sorted_bed],
               stdout=open(merged_bed, 'w'), check=True)

def n(s):
    return int(subprocess.check_output(s, shell=True).decode().strip() or '0')

tr = n(f"awk '{{sum+=$3-$2}} END {{print sum}}' {TRUTH_BED}")
md = n(f"awk '{{sum+=$3-$2}} END {{print sum}}' {merged_bed}")
tp = n(f"{BEDTOOLS} intersect -a {TRUTH_BED} -b {merged_bed} | "
       f"awk '{{sum+=$3-$2}} END {{print sum}}'")

if md == 0 or tr == 0:
    print(f"[{LABEL}] FAIL: empty mask or truth")
    sys.exit(1)

r = tp / tr; p = tp / md; f = 2 * r * p / (r + p) if (r + p) > 0 else 0
print(f"[{LABEL}] recall={r:.4f} precision={p:.4f} F1={f:.4f} mask={md/1e6:.2f}Mb tp={tp/1e6:.2f}Mb")
