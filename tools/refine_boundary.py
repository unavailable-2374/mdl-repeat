#!/usr/bin/env python3
"""Boundary-only consensus refinement for mdl-repeat output.

Iterative post-processing that extends each family's consensus boundaries
based on rmblastn-recruited genomic instances. NEVER modifies the internal
consensus — preserves RepeatMasker matching specificity.

Validated 2026-05-02 on TAIR10 nuclear:
  + 7 iterations → bp-F1 +2.00 pp (0.5793 → 0.5993)
  + bp-recall +2.54 pp; bp-precision +0.17 pp (extensions are real TE bp)

Algorithm
---------
For each family with ≥ MIN_INSTANCES_FOR_EXT rmblastn hits at ≥ IDENTITY_THRESHOLD:
  5' side:
    1. Find hits anchored at 5' end (qstart ≤ ANCHOR_TOLERANCE)
    2. Extract EXT_LOOKAHEAD bp upstream of each hit's sstart on the genome
       (strand-aware: for - strand, take downstream of send and reverse-complement)
    3. Going LEFT from consensus position 0:
       - At each offset, check if ≥ MIN_OCCUPANCY of recruited flanks have
         a non-N base; if not → stop
       - Compute per-position majority base; if majority ≥ 50% → append
       - Cap at MAX_EXT_PER_SIDE
  3' side: symmetric (qend ≥ qlen - ANCHOR_TOLERANCE + 1, downstream of send)
  Output: (5' extension) + (original consensus) + (3' extension)

Iterate by feeding output back as input — each iteration extends further as
new copies become "anchored" after previous extensions.

Why boundary-only (not whole-consensus)
---------------------------------------
We tested 6 variants of whole-consensus refinement (hmmemit, MAFFT majority,
augmented library) — all caused 15+ pp bp-F1 drops via RepeatMasker. Diagnosis:
mdl-repeat's column-majority consensus is well-tuned for cross-match anchoring;
replacing it with "averaged" sequences breaks RM extension specificity. ONLY
boundary-only extension preserves the central consensus byte-for-byte AND
adds value.

Why iterate
-----------
Per-iteration TAIR10 gains: 1.06 → 0.35 → 0.22 → 0.13 → 0.09 → 0.08 → 0.07 (asymp).
Each iter exposes new copies that anchor at the new (extended) boundary.

Dependencies (in EDTA conda env or equivalent)
----------------------------------------------
  - rmblastn (RepeatMasker BLAST)
  - makeblastdb
  - samtools

Usage
-----
  python3 tools/refine_boundary.py <library.fa> <genome.fa> <output_dir> [N_iters]
  → writes refined_v6.fa, log_v6.tsv into output_dir per iteration
  → if N_iters > 1: outputs refined_iter{i}.fa for each iteration

Defaults work for plant genomes; tune via env vars if needed:
  REFINE_IDENTITY=70  REFINE_LOOKAHEAD=100  REFINE_MAX_EXT=100
"""
import os, sys, subprocess, tempfile, shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict, Counter

if len(sys.argv) < 4:
    print(__doc__, file=sys.stderr)
    sys.exit(1)

INPUT_FA = sys.argv[1]
GENOME_FA = sys.argv[2]
OUTDIR = sys.argv[3]
N_ITERS = int(sys.argv[4]) if len(sys.argv) > 4 else 1
os.makedirs(OUTDIR, exist_ok=True)

# Tools — locate from PATH or env override
RMBLASTN = os.environ.get('RMBLASTN', 'rmblastn')
MAKEBLASTDB = os.environ.get('MAKEBLASTDB', 'makeblastdb')
SAMTOOLS = os.environ.get('SAMTOOLS', 'samtools')

# Tunable parameters
IDENTITY_THRESHOLD = float(os.environ.get('REFINE_IDENTITY', '70'))
COV_THRESHOLD = float(os.environ.get('REFINE_COVERAGE', '0.40'))
MIN_INSTANCES_FOR_EXT = int(os.environ.get('REFINE_MIN_INSTANCES', '3'))
ANCHOR_TOLERANCE = int(os.environ.get('REFINE_ANCHOR_TOL', '5'))
EXT_LOOKAHEAD = int(os.environ.get('REFINE_LOOKAHEAD', '100'))
MIN_OCCUPANCY = float(os.environ.get('REFINE_MIN_OCC', '0.50'))
MAX_EXT_PER_SIDE = int(os.environ.get('REFINE_MAX_EXT', '100'))
N_THREADS = int(os.environ.get('REFINE_THREADS', '4'))

COMP = {'A':'T','T':'A','C':'G','G':'C','N':'N','a':'t','t':'a','c':'g','g':'c','n':'n'}
def revcomp(s):
    return ''.join(COMP.get(c, 'N') for c in reversed(s))

def parse_fasta(path):
    seqs = {}
    cur = None; buf = []
    for ln in open(path):
        if ln.startswith('>'):
            if cur: seqs[cur] = ''.join(buf)
            cur = ln[1:].rstrip().split()[0]; buf = []
        else: buf.append(ln.rstrip())
    if cur: seqs[cur] = ''.join(buf)
    return seqs

def extend_one_side(qid, fam_hits, side, qlen, genome_fa):
    if side == '5':
        anchored = [h for h in fam_hits if h['qstart'] <= ANCHOR_TOLERANCE]
    else:
        anchored = [h for h in fam_hits if h['qend'] >= qlen - ANCHOR_TOLERANCE + 1]
    if len(anchored) < MIN_INSTANCES_FOR_EXT:
        return ''

    flanks = []
    for h in anchored[:30]:
        sid, sstart, send, strand = h['sid'], h['sstart'], h['send'], h['strand']
        if side == '5':
            if strand == '+':
                ext_start = max(1, sstart - EXT_LOOKAHEAD)
                ext_end = sstart - 1
                rc = False
            else:
                ext_start = send + 1
                ext_end = send + EXT_LOOKAHEAD
                rc = True
        else:
            if strand == '+':
                ext_start = send + 1
                ext_end = send + EXT_LOOKAHEAD
                rc = False
            else:
                ext_start = max(1, sstart - EXT_LOOKAHEAD)
                ext_end = sstart - 1
                rc = True
        if ext_end < ext_start: continue
        try:
            out = subprocess.check_output(
                [SAMTOOLS, 'faidx', genome_fa, f"{sid}:{ext_start}-{ext_end}"],
                timeout=5, stderr=subprocess.DEVNULL).decode()
            seq = ''.join(out.strip().split('\n')[1:]).upper()
            if rc: seq = revcomp(seq)
            if seq: flanks.append(seq)
        except Exception:
            continue

    if len(flanks) < MIN_INSTANCES_FOR_EXT:
        return ''

    extension = []
    for offset in range(MAX_EXT_PER_SIDE):
        bases = []
        for f in flanks:
            idx = (len(f) - 1 - offset) if side == '5' else offset
            if 0 <= idx < len(f):
                bases.append(f[idx])
        if len(bases) < MIN_INSTANCES_FOR_EXT:
            break
        non_n = [b for b in bases if b != 'N']
        if len(non_n) / len(flanks) < MIN_OCCUPANCY:
            break
        c = Counter(non_n)
        if not c: break
        best, cnt = c.most_common(1)[0]
        if cnt / len(non_n) < 0.5:
            break
        extension.append(best)
    return ''.join(reversed(extension)) if side == '5' else ''.join(extension)

def refine_one(args):
    qid, fam_hits, original, genome_fa = args
    if len(fam_hits) < MIN_INSTANCES_FOR_EXT:
        return qid, original, 0, 0, 'too_few_total'
    qlen = fam_hits[0]['qlen']
    ext5 = extend_one_side(qid, fam_hits, '5', qlen, genome_fa)
    ext3 = extend_one_side(qid, fam_hits, '3', qlen, genome_fa)
    new_seq = ext5 + original + ext3
    return qid, new_seq, len(ext5), len(ext3), 'extended' if (ext5 or ext3) else 'no_extension'

def run_one_iteration(input_fa, genome_fa, outdir, iter_label=''):
    families = parse_fasta(input_fa)
    print(f"[{iter_label}] Loaded {len(families)} families from {input_fa}", file=sys.stderr)

    db_path = os.path.join(outdir, 'genome_db')
    if not os.path.exists(db_path + '.nhr'):
        subprocess.run([MAKEBLASTDB, '-in', genome_fa, '-dbtype', 'nucl',
                        '-out', db_path, '-parse_seqids'],
                       check=True, capture_output=True)

    blast_out = os.path.join(outdir, f'lib_vs_genome{iter_label}.tsv')
    if not os.path.exists(blast_out):
        print(f"[{iter_label}] rmblastn ≥{IDENTITY_THRESHOLD}% ...", file=sys.stderr)
        subprocess.run([
            RMBLASTN, '-query', input_fa, '-db', db_path,
            '-outfmt', '6 qseqid sseqid pident length qlen slen sstart send qstart qend evalue',
            '-evalue', '1e-10', '-num_threads', str(N_THREADS),
            '-perc_identity', str(IDENTITY_THRESHOLD),
            '-gapopen', '5', '-gapextend', '2',
            '-dust', 'no', '-soft_masking', 'false',
            '-out', blast_out,
        ], check=True, capture_output=True)
    print(f"[{iter_label}] rmblastn done: {sum(1 for _ in open(blast_out))} hits", file=sys.stderr)

    hits = defaultdict(list)
    with open(blast_out) as f:
        for ln in f:
            p = ln.rstrip().split('\t')
            if len(p) < 11: continue
            qid, sid = p[0], p[1]
            pident = float(p[2])
            alen, qlen, slen = int(p[3]), int(p[4]), int(p[5])
            sstart, send = int(p[6]), int(p[7])
            qstart, qend = int(p[8]), int(p[9])
            if alen / qlen < COV_THRESHOLD: continue
            strand = '+' if sstart < send else '-'
            if strand == '-':
                sstart, send = send, sstart
            hits[qid].append({
                'sid': sid, 'pident': pident,
                'sstart': sstart, 'send': send,
                'qstart': qstart, 'qend': qend,
                'qlen': qlen, 'slen': slen,
                'strand': strand, 'alen': alen,
            })

    work = [(qid, hits.get(qid, []), seq, genome_fa) for qid, seq in families.items()]
    results = {}
    log = []
    done = 0
    with ProcessPoolExecutor(max_workers=N_THREADS) as exe:
        futures = {exe.submit(refine_one, w): w[0] for w in work}
        for fut in as_completed(futures):
            try:
                qid, new_seq, e5, e3, status = fut.result(timeout=120)
                results[qid] = new_seq
                log.append((qid, status, e5, e3, len(families[qid]), len(new_seq)))
            except Exception as e:
                qid = futures[fut]
                results[qid] = families[qid]
                log.append((qid, 'error', 0, 0, len(families[qid]), len(families[qid])))
            done += 1
            if done % 1000 == 0:
                print(f"[{iter_label}]   done {done}/{len(work)}", file=sys.stderr)

    out_fa = os.path.join(outdir, f'refined{iter_label}.fa')
    with open(out_fa, 'w') as fh:
        for qid in families:
            fh.write(f">{qid}\n{results.get(qid, families[qid])}\n")
    log_path = os.path.join(outdir, f'log{iter_label}.tsv')
    with open(log_path, 'w') as fh:
        fh.write('family_id\tstatus\text5\text3\torig_len\tnew_len\n')
        for row in log: fh.write('\t'.join(str(x) for x in row) + '\n')

    n_ext = sum(1 for _, status, _, _, _, _ in log if status == 'extended')
    total5 = sum(e5 for _, _, e5, _, _, _ in log)
    total3 = sum(e3 for _, _, _, e3, _, _ in log)
    print(f"[{iter_label}] {n_ext} families extended; +{total5} 5' bp, +{total3} 3' bp", file=sys.stderr)
    return out_fa

cur_input = INPUT_FA
for it in range(1, N_ITERS + 1):
    label = f'_iter{it}' if N_ITERS > 1 else ''
    cur_input = run_one_iteration(cur_input, GENOME_FA, OUTDIR, label)
print(f"\nFinal: {cur_input}", file=sys.stderr)
